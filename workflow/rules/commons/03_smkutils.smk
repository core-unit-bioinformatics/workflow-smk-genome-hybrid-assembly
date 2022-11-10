
localrules:
    dump_config,
    create_manifest,


rule accounting_file_md5_size:
    """
    Compute MD5 checksum and file size
    for accounted files (inputs, references, results).
    For convenience, this rule also computes the file
    size to avoid this overhead at some other place
    of the pipeline.
    """
    input:
        source=load_file_by_path_id,
    output:
        md5=DIR_PROC.joinpath(
            ".accounting", "checksums", "{file_type}", "{file_name}.{path_id}.md5"
        ),
        file_size=DIR_PROC.joinpath(
            ".accounting", "file_sizes", "{file_type}", "{file_name}.{path_id}.bytes"
        ),
    benchmark:
        DIR_RSRC.joinpath(
            ".accounting", "checksums", "{file_type}", "{file_name}.{path_id}.md5.rsrc"
        )
    wildcard_constraints:
        file_type="(" + "|".join(sorted(ACCOUNTING_FILES.keys())) + ")",
    resources:
        time_hrs=lambda wildcards, attempt: 1 * attempt,
        mem_gb=lambda wildcards, attempt: 1 * attempt,
    shell:
        "md5sum {input.source} > {output.md5}"
        " && "
        "stat -c %s {input.source} > {output.file_size}"


rule accounting_file_sha256:
    """
    Compute SHA256 checksum (same as for MD5)
    """
    input:
        source=load_file_by_path_id,
    output:
        sha256=DIR_PROC.joinpath(
            ".accounting", "checksums", "{file_type}", "{file_name}.{path_id}.sha256"
        ),
    benchmark:
        DIR_RSRC.joinpath(
            ".accounting",
            "checksums",
            "{file_type}",
            "{file_name}.{path_id}.sha256.rsrc",
        )
    wildcard_constraints:
        file_type="(" + "|".join(sorted(ACCOUNTING_FILES.keys())) + ")",
    resources:
        time_hrs=lambda wildcards, attempt: 1 * attempt,
        mem_gb=lambda wildcards, attempt: 1 * attempt,
    shell:
        "sha256sum {input.source} > {output.sha256}"


rule dump_config:
    output:
        RUN_CONFIG_RELPATH,
    params:
        acc_in=lambda wildcards, output: register_input(output, allow_non_existing=True),
    run:
        import yaml

        runinfo = {"_timestamp": get_timestamp(), "_username": get_username()}

        try:
            git_labels = collect_git_labels()
        except Exception as err:
            logerr(f"dump_config: collect_git_labels() raised exception: {err}")
            raise
        for label, value in git_labels:
            runinfo[f"_{label}"] = value
        # add complete Snakemake config
        runinfo.update(config)
        for special_key in ["devmode", "resetacc"]:
            try:
                del runinfo[special_key]
            except KeyError:
                pass

        with open(RUN_CONFIG_RELPATH, "w", encoding="ascii") as cfg_dump:
            yaml.dump(runinfo, cfg_dump, allow_unicode=False, encoding="ascii")
        # END OF RUN BLOCK



rule create_manifest:
    input:
        manifest_files=load_accounting_information,
    output:
        manifest=MANIFEST_RELPATH,
    run:
        import fileinput
        import collections
        import pandas

        process_accounting_files = {}
        for accounting_file, file_path in ACCOUNTING_FILES.items():
            if not file_path.is_file():
                if VERBOSE:
                    warn_msg = f"Warning: accounting file of type '{account_file}' not in use."
                    logerr(warn_msg)
                continue
            process_accounting_files[accounting_file] = file_path

        if len(process_accounting_files) == 0:
            err_msg = "No accounting files marked are in use.\n"
            err_msg += "This means one of two things:\n"
            err_msg += "1) Your workflow does not consume input, does not use "
            err_msg += "any reference information and also does not produce output.\n"
            err_msg += "Really? Are you sure?\n"
            err_msg += "2) You did not annotate the workflow rules with:\n"
            err_msg += "commons/02_pyutils.smk::register_input()\n"
            err_msg += "commons/02_pyutils.smk::register_result()\n"
            err_msg += "commons/02_pyutils.smk::register_reference()\n"
            err_msg += "Please rerun the workflow twice in dry run mode...\n"
            err_msg += "snakemake --dry-run (or: -n) [...other options...]\n"
            err_msg += "...after fixing that.\n"
            err_msg += "However, if you are sure that is correct, remove the entry\n"
            err_msg += "MANIFEST_RELPATH from the rule 'run_all' in the 'workflow/Snakefile'\n"
            err_msg += "and restart the workflow."
            logerr(err_msg)
            raise RuntimeError(
                "No accounts: cannot proceed with workflow execution w/o accouting files."
            )

        if len(input.manifest_files) == 0:
            # after checking accounting files in use, this
            # should rather point to an error
            warn_msg = "No files recorded for inclusion in workflow manifest.\n"
            warn_msg += "Are you sure you did not forget annotating rules with:\n"
            warn_msg += "commons/02_pyutils.smk::register_input()\n"
            warn_msg += "commons/02_pyutils.smk::register_result()\n"
            logerr(warn_msg)

        records = collections.defaultdict(dict)
        for line in fileinput.input(process_accounting_files.values(), mode="r"):
            path_id, path_record = process_accounting_record(line)
            records[path_id].update(path_record)

        df = pandas.DataFrame.from_records(list(records.values()))
        if df.empty:
            logerr("Manifest DataFrame is empty - aborting")
            raise RuntimeError("Manifest DataFrame is empty")

        df.sort_values(["file_category", "file_name"], ascending=True, inplace=True)
        reordered_columns = [
            "file_name",
            "file_category",
            "file_size_bytes",
            "file_checksum_md5",
            "file_checksum_sha256",
            "file_path",
            "path_id",
        ]
        assert all(c in df.columns for c in reordered_columns)
        df = df[reordered_columns]
        df.to_csv(output.manifest, header=True, index=False, sep="\t")
        # END OF RUN BLOCK
