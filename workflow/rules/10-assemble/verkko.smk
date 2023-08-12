localrules: verkko_trio_samples
rule verkko_trio_samples:
    """
    This rule sets the Verkko option
    --lsf
    to force-run Verkko in some cluster
    mode. It is vital for proper execution
    that ALL Snakemake profile parameters
    are properly set in the profile supplied
    via --snakeopts to override the defaults.
    """
    input:
        hifi = lambda wildcards: MAP_SAMPLE_TO_INPUT_FILES[wildcards.sample]["hifi"],
        nano = lambda wildcards: MAP_SAMPLE_TO_INPUT_FILES[wildcards.sample]["ont"],
        hap1_db = lambda wildcards: MAP_SAMPLE_TO_INPUT_FILES[wildcards.sample]["hap1"],
        hap2_db = lambda wildcards: MAP_SAMPLE_TO_INPUT_FILES[wildcards.sample]["hap2"],
        profile = ancient(config["verkko_smk_profile"]),
    output:
        done = DIR_PROC.joinpath("10-assemble/verkko/{sample}.ps-trio.ok")
    log:
        DIR_LOG.joinpath("10-assemble/verkko/{sample}.ps-trio.log")
    benchmark:
        DIR_RSRC.joinpath("10-assemble/verkko/{sample}.ps-trio.rsrc")
    conda:
        DIR_ENVS.joinpath("verkko.yaml")
    wildcard_constraints:
        sample = CONSTRAINT_TRIO_SAMPLES
    resources:
        mbg_rsrc=lambda wildcards, attempt: increase_mbg_resources(attempt),
        pop_rsrc=lambda wildcards, attempt: increase_process_ont_resources(attempt),
        lay_rsrc=lambda wildcards, attempt: increase_layout_contigs_resources(attempt),
    params:
        dryrun="--dry-run" if VERKKO_DRY_RUN else "",
        check=lambda wildcards, output: "" if VERKKO_DRY_RUN else f" && touch {output.done}",
        acc_in=lambda wildcards, input: register_input(input.hifi, input.nano),
        screen=lambda wildcards: assemble_verkko_screen_string(),
        wd=lambda wildcards, output: pathlib.Path(output.done).with_suffix(".wd"),
    shell:
        "/usr/bin/time -v "
        "verkko --lsf "
            "--python `which python` "
            "--mbg `which MBG` "
            "--graphaligner `which GraphAligner` "
            "--mashmap `which mashmap` "
            "--hifi {input.hifi} "
            "--nano {input.nano} "
            "{params.screen} "
            "-d {params.wd} "
            "{resources.mbg_rsrc} "
            "{resources.pop_rsrc} "
            "{resources.lay_rsrc} "
            "--hap-kmers {input.hap1_db} {input.hap2_db} trio "
            "--snakeopts \"--profile $PWD/{input.profile} {params.dryrun}\" "
            "&> {log} {params.check}"


localrules: verkko_unphased_samples
rule verkko_unphased_samples:
    """
    This rule sets the Verkko option
    --lsf
    to force-run Verkko in some cluster
    mode. It is vital for proper execution
    that ALL Snakemake profile parameters
    are properly set in the profile supplied
    via --snakeopts to override the defaults.
    """
    input:
        hifi = lambda wildcards: MAP_SAMPLE_TO_INPUT_FILES[wildcards.sample]["hifi"],
        nano = lambda wildcards: MAP_SAMPLE_TO_INPUT_FILES[wildcards.sample]["ont"],
        profile = ancient(config["verkko_smk_profile"]),
    output:
        done = DIR_PROC.joinpath("assemblies/verkko/{sample}.ps-none.ok")
    log:
        DIR_LOG.joinpath("assemblies/verkko/{sample}.ps-none.log")
    benchmark:
        DIR_RSRC.joinpath("assemblies/verkko/{sample}.ps-none.rsrc")
    conda:
        DIR_ENVS.joinpath("verkko.yaml")
    wildcard_constraints:
        sample = CONSTRAINT_UNPHASED_SAMPLES
    resources:
        mbg_rsrc=lambda wildcards, attempt: increase_mbg_resources(attempt),
        pop_rsrc=lambda wildcards, attempt: increase_process_ont_resources(attempt),
        lay_rsrc=lambda wildcards, attempt: increase_layout_contigs_resources(attempt),
    params:
        dryrun = "--dry-run" if VERKKO_DRY_RUN else "",
        check = lambda wildcards, output: "" if VERKKO_DRY_RUN else f" && touch {output.done}",
        acc_in=lambda wildcards, input: register_input(input.hifi, input.nano),
        screen=lambda wildcards: assemble_verkko_screen_string(),
        wd=lambda wildcards, output: pathlib.Path(output.done).with_suffix(".wd"),
    shell:
        "/usr/bin/time -v "
        "verkko --lsf "
            "--python `which python` "
            "--mbg `which MBG` "
            "--graphaligner `which GraphAligner` "
            "--mashmap `which mashmap` "
            "--hifi {input.hifi} "
            "--nano {input.nano} "
            "{params.screen} "
            "-d {params.wd} "
            "{resources.mbg_rsrc} "
            "{resources.pop_rsrc} "
            "{resources.lay_rsrc} "
            "--snakeopts \"--profile $PWD/{input.profile} {params.dryrun}\" "
            "&> {log} {params.check}"


localrules: verkko_strandseq_samples
rule verkko_strandseq_samples:
    """
    This rule sets the Verkko option
    --lsf
    to force-run Verkko in some cluster
    mode. It is vital for proper execution
    that ALL Snakemake profile parameters
    are properly set in the profile supplied
    via --snakeopts to override the defaults.
    """
    input:
        unphased = DIR_PROC.joinpath("assemblies/verkko/{sample}.ps-none.ok"),
        hifi = lambda wildcards: MAP_SAMPLE_TO_INPUT_FILES[wildcards.sample]["hifi"],
        nano = lambda wildcards: MAP_SAMPLE_TO_INPUT_FILES[wildcards.sample]["ont"],
        paths = lambda wildcards: MAP_SAMPLE_TO_INPUT_FILES[wildcards.sample]["phasing_paths"],
        profile = ancient(config["verkko_smk_profile"]),
    output:
        done = DIR_PROC.joinpath("10-assemble/verkko/{sample}.ps-sseq.wait")
    log:
        DIR_LOG.joinpath("10-assemble/verkko/{sample}.ps-sseq.log")
    benchmark:
        DIR_RSRC.joinpath("10-assemble/verkko/{sample}.ps-sseq.rsrc")
    conda:
        DIR_ENVS.joinpath("verkko.yaml")
    wildcard_constraints:
        sample = CONSTRAINT_SSEQ_SAMPLES
    resources:
        mbg_rsrc=lambda wildcards, attempt: increase_mbg_resources(attempt),
        pop_rsrc=lambda wildcards, attempt: increase_process_ont_resources(attempt),
        lay_rsrc=lambda wildcards, attempt: increase_layout_contigs_resources(attempt),
    params:
        dryrun = "--dry-run" if VERKKO_DRY_RUN else "",
        check = lambda wildcards, output: "" if VERKKO_DRY_RUN else f" && touch {output.done}",
        acc_in=lambda wildcards, input: register_input(input.hifi, input.nano),
        screen=lambda wildcards: assemble_verkko_screen_string(),
        assm_dir=lambda wildcards, input: pathlib.Path(input.unphased).with_suffix(".wd"),
        phasing_dir=lambda wildcards, output: pathlib.Path(output.done).with_suffix(".wd")
    shell:
        "/usr/bin/time -v "
        "verkko --lsf "
            "--python `which python` "
            "--mbg `which MBG` "
            "--graphaligner `which GraphAligner` "
            "--mashmap `which mashmap` "
            "--hifi {input.hifi} "
            "--nano {input.nano} "
            "{params.screen} "
            "-d {params.phasing_dir} "
            "--assembly {params.assm_dir} "
            "--paths {input.paths} "
            "{resources.mbg_rsrc} "
            "{resources.pop_rsrc} "
            "{resources.lay_rsrc} "
            "--snakeopts \"--profile $PWD/{input.profile} {params.dryrun}\" "
            "&> {log} {params.check}"


rule confirm_path_consistency:
    """This rule is a (temporary?)
    sanity check to confirm that the
    phasing paths generated with the Strand-seq
    phasing pipeline (i.e., external input)
    are identical to the ones then used
    for restarting the Verkko assembly
    process with the --paths argument.
    """
    input:
        check_file = rules.verkko_strandseq_samples.output.done,
        original = lambda wildcards: MAP_SAMPLE_TO_INPUT_FILES[wildcards.sample]["phasing_paths"],
    output:
        confirmed = DIR_PROC.joinpath(
            "10-assemble", "verkko", "{sample}.ps-sseq.paths.ok"
        )
    run:
        import hashlib as hl
        import pathlib as pl

        copied_paths = pl.Path(input.check_file).with_suffix(".wd").joinpath(
            "6-layoutContigs", "consensus_paths.txt"
        )
        assert copied_paths.is_file(), f"Verkko run incomplete: {wildcards.sample}.ps-sseq"

        original_md5 = hl.md5(open(input.original, "rb").read()).hexdigest()
        copied_md5 = hl.md5(open(copied_paths, "rb").read()).hexdigest()

        if original_md5 != copied_md5:
            logerr(f"Verkko phasing paths inconsistent: {input.original} vs {copied_paths}")
            raise RuntimeError("10-assemble::verkko.smk::confirm_path_consistency - no MD5 match")

        with open(output.confirmed, "w") as dump:
            _ = dump.write(
                f"{get_timestamp()}\n{input.original}\n{copied_paths}\n{original_md5}\n"
            )
    # END OF RUN BLOCK


rule run_verkko_unphased_samples:
    input:
        assemblies = expand(
            DIR_PROC.joinpath("assemblies/verkko/{sample}.ps-none.ok"),
            sample=UNPHASED_SAMPLES
        ),


rule run_verkko_trio_samples:
    input:
        assemblies = expand(
            DIR_PROC.joinpath("10-assemble/verkko/{sample}.ps-trio.ok"),
            sample=TRIO_SAMPLES
        ),


rule run_verkko_sseq_samples:
    input:
        assemblies = expand(
            DIR_PROC.joinpath("10-assemble/verkko/{sample}.ps-sseq.ok"),
            sample=SSEQ_SAMPLES
        ),
