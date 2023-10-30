
rule split_verkko_posthoc_phased_fasta:
    """
    This is a temporary rule to deal with
    the issue that Verkko does not split
    the main fasta into hap1/hap2 and
    unassigned if Verkko is restarted
    with the --paths argument.
    This should be fixed in some future
    release of Verkko.

    See gh#170
    """
    input:
        wait_file = DIR_PROC.joinpath("10-assemble/verkko/{sample}.{phasing_state}.wait"),
        path_check = DIR_PROC.joinpath("10-assemble/verkko/{sample}.{phasing_state}.paths.ok")
    output:
        check_file = DIR_PROC.joinpath("10-assemble/verkko/{sample}.{phasing_state}.ok")
    wildcard_constraints:
        phasing_state = "ps-sseq"
    resources:
        mem_mb=lambda wildcards, attempt: 8192 * attempt
    run:
        import io
        import os
        import pathlib as pl

        def is_empty(file_path):
            return os.stat(file_path).st_size == 0

        def dump_fasta_partition(part_name, part_seq_buffer, verkko_wd):
            part_file = verkko_wd.joinpath(f"assembly.{part_name}.fasta")
            assert part_file.is_file()
            assert is_empty(part_file)
            with open(part_file, "w") as fasta:
                _ = fasta.write(part_seq_buffer.getvalue())
            return

        verkko_run_wd = pl.Path(input.wait_file).with_suffix(".wd")
        assert verkko_run_wd.is_dir()

        main_assembly_file = verkko_run_wd.joinpath("assembly.fasta")
        assert main_assembly_file.is_file()
        assert not is_empty(main_assembly_file)

        expected_partitions = ["haplotype1", "haplotype2", "unassigned"]
        processed_partitions = set()
        active_partition = None
        partition_buffer = io.StringIO()
        num_contigs = 0
        with open(main_assembly_file, "r") as fasta:
            for line in fasta:
                if line.startswith(">"):
                    contig_partition = line.strip()[1:].split("-")[0]
                    num_contigs += 1
                    assert contig_partition in expected_partitions
                    if processed_partitions and contig_partition not in processed_partitions:
                        dump_fasta_partition(active_partition, partition_buffer, verkko_run_wd)
                        processed_partitions.add(active_partition)
                        partition_buffer = io.StringIO()
                        active_partition = contig_partition
                    active_partition = contig_partition
                    processed_partitions.add(active_partition)
                    partition_buffer.write(line)
                else:
                    partition_buffer.write(line)

        dump_fasta_partition(active_partition, partition_buffer, verkko_run_wd)
        processed_partitions.add(active_partition)
        assert sorted(processed_partitions) == expected_partitions, (
            f"Partitions: processed {processed_partitions} - expected {expected_partitions}"
        )

        with open(output.check_file, "w") as dump:
            _ = dump.write(f"Processed contigs: {num_contigs}\n")
    # END OF RUN BLOCK


localrules: collect_verkko_output_files
rule collect_verkko_output_files:
    input:
        check_file = DIR_PROC.joinpath("10-assemble/verkko/{sample}.{phasing_state}.ok")
    output:
        file_collection = DIR_PROC.joinpath("10-assemble/verkko/{sample}.{phasing_state}.output.json")
    conda:
        DIR_ENVS.joinpath("pygraph.yaml")
    params:
        script = find_script("collect_verkko_output"),
        verkko_wd = lambda wildcards, input: pathlib.Path(input.check_file).with_suffix(".wd"),
        base_dir = WORKDIR
    shell:
        "{params.script} --verkko-wd {params.verkko_wd} --base-dir {params.base_dir} "
            "--delete-logs --sample {wildcards.sample} --output {output.file_collection}"


################################
### 2023-08-03
# Verkko assembly output consists
# of duplicated sequences at least
# between disconnected and rDNA and
# disconnected and EBV. Unclear if
# discarding that would result in
# information loss, but for the
# time being, it is assumed that
# this is not the case. New rules
# added to deduplicate the output
# at least in the non-primary
# output files ("scraps").
###
#################################

rule compute_verkko_assembled_sequence_id:
    """
    The output computed here will eventually
    be discarded, and is only used to identify
    identical sequences among all Verkko output
    files (by sequence hash).
    """
    input:
        file_collection = DIR_PROC.joinpath("10-assemble/verkko/{sample}.{phasing_state}.output.json")
    output:
        stats = DIR_PROC.joinpath(
            "30-postprocess", "deduplicate", "seq_ids",
            "{sample}.{phasing_state}.verkko-asm-{asmtype}.statistics.tsv.gz"
        ),
    wildcard_constraints:
        asmtype = "(hap1|hap2|unassigned|disconnected|rdna|ebv|mito|wg)"
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    threads: CPU_LOW
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt,
        time_hrs=lambda wildcards, attempt: attempt,
    params:
        script=find_script("seqstats"),
        assembly=lambda wildcards, input: get_verkko_output(
            input.file_collection, f"{wildcards.asmtype}_fasta"
        ),
    shell:
        "{params.script} --cores {threads} "
        "--no-homopolymer-runs "
        "--no-short-tandem-repeats "
        "--no-sequence-composition "
        "--output-statistics {output.stats} "
        "--input-files {params.assembly}"


rule prepare_verkko_dup_seq_report:
    input:
        seq_stats = lambda wildcards: expand(
            rules.compute_verkko_assembled_sequence_id.output.stats,
            asmtype=get_verkko_asm_units(wildcards.phasing_state),
            allow_missing=True
        )
    output:
        report = DIR_RES.joinpath(
            "reports", "seq_dedup", "{sample}.{phasing_state}.verkko-dup-seq.tsv.gz"
        ),
        summary = DIR_RES.joinpath(
            "reports", "seq_dedup", "{sample}.{phasing_state}.verkko-dup-seq.summary.tsv"
        ),
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        script = find_script("asm_dupseq_report")
    shell:
        "{params.script} --input {input.seq_stats} --report {output.report} "
            "--summary {output.summary}"


rule filter_verkko_dup_sequences:
    input:
        report = rules.prepare_verkko_dup_seq_report.output.report,
        file_collection = rules.collect_verkko_output_files.output.file_collection
    output:
        asm_unit = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "{sample}.{phasing_state}.asm-{asm_unit}.fasta.gz"
        ),
        fai = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "{sample}.{phasing_state}.asm-{asm_unit}.fasta.gz.fai"
        ),
        gzi = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "{sample}.{phasing_state}.asm-{asm_unit}.fasta.gz.gzi"
        ),
    log:
        DIR_LOG.joinpath(
            "30-postprocess", "verkko",
            "{sample}.{phasing_state}.asm-{asm_unit}.dedup.log"
        ),
    wildcard_constraints:
        asm_unit = "(hap1|hap2|unassigned|disconnected|rdna|ebv|mito|wg)"
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        script = find_script("filter_verkko_dupseq"),
        fasta = lambda wildcards, input: get_verkko_output(
            input.file_collection, f"{wildcards.asm_unit}_fasta"
        ),
        acc_res=lambda wildcards, output: register_result(output)
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        "rm -f {output.asm_unit}.EMPTY ; "
        "if [ -s {params.fasta} ] ; then "
        "{{ "
        "{params.script} --input {params.fasta} --report {input.report} --verbose 2> {log}"
            " | "
        "bgzip > {output.asm_unit}"
            " && "
        "samtools faidx {output.asm_unit} ; "
        "}} else {{ "
        "touch {output.asm_unit} && touch {output.fai} "
        "&& touch {output.gzi} && touch {output.asm_unit}.EMPTY ; "
        "}} fi ;"


localrules: copy_verkko_exemplar_sequences
rule copy_verkko_exemplar_sequences:
    input:
        file_collection = rules.collect_verkko_output_files.output.file_collection
    output:
        ex_seq = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "{sample}.{phasing_state}.exemplar-{asm_unit}.fasta.gz"
        ),
        fai = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "{sample}.{phasing_state}.exemplar-{asm_unit}.fasta.gz.fai"
        ),
        gzi = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "{sample}.{phasing_state}.exemplar-{asm_unit}.fasta.gz.gzi"
        ),
    wildcard_constraints:
        asm_unit = "(mito|ebv|rdna)"
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        fasta = lambda wildcards, input: get_verkko_output(
            input.file_collection, f"{wildcards.asm_unit}_repr"
        ),
        acc_res=lambda wildcards, output: register_result(output)
    shell:
        "rm -f {output.ex_seq}.EMPTY ; "
        "if [ -s {params.fasta} ] ; then "
        "{{ "
        "cat {params.fasta} | bgzip > {output.ex_seq}"
            " && "
        "samtools faidx {output.ex_seq} ; "
        "}} else {{ "
        "touch {output.ex_seq} && touch {output.fai} "
        "&& touch {output.gzi} touch {output.ex_seq}.EMPTY ; "
        "}} fi;"


rule postprocess_verkko_unphased_samples:
    input:
        exemplars = expand(
            rules.copy_verkko_exemplar_sequences.output.ex_seq,
            phasing_state=["ps-none"],
            asm_unit=["rdna", "ebv", "mito"],
            sample=SAMPLES
        ),
        asm_units = expand(
            rules.filter_verkko_dup_sequences.output.asm_unit,
            phasing_state=["ps-none"],
            asm_unit=get_verkko_asm_units("ps-none"),
            sample=SAMPLES
        )


rule postprocess_verkko_sseq_samples:
    input:
        exemplars = expand(
            rules.copy_verkko_exemplar_sequences.output.ex_seq,
            phasing_state=["ps-sseq"],
            asm_unit=["rdna", "ebv", "mito"],
            sample=SSEQ_SAMPLES
        ),
        asm_units = expand(
            rules.filter_verkko_dup_sequences.output.asm_unit,
            phasing_state=["ps-sseq"],
            asm_unit=get_verkko_asm_units("ps-sseq"),
            sample=SSEQ_SAMPLES
        )


rule postprocess_verkko_trio_samples:
    input:
        exemplars = expand(
            rules.copy_verkko_exemplar_sequences.output.ex_seq,
            phasing_state=["ps-trio"],
            asm_unit=["rdna", "ebv", "mito"],
            sample=TRIO_SAMPLES
        ),
        asm_units = expand(
            rules.filter_verkko_dup_sequences.output.asm_unit,
            phasing_state=["ps-trio"],
            asm_unit=get_verkko_asm_units("ps-trio"),
            sample=TRIO_SAMPLES
        )


rule postprocess_verkko_hic_samples:
    input:
        exemplars = expand(
            rules.copy_verkko_exemplar_sequences.output.ex_seq,
            phasing_state=["ps-hic"],
            asm_unit=["rdna", "ebv", "mito"],
            sample=HIC_SAMPLES
        ),
        asm_units = expand(
            rules.filter_verkko_dup_sequences.output.asm_unit,
            phasing_state=["ps-hic"],
            asm_unit=get_verkko_asm_units("ps-hic"),
            sample=HIC_SAMPLES
        )
