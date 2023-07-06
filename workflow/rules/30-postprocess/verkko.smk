
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
        wait_file = DIR_PROC.joinpath("assemblies/verkko/{sample}.{phasing_state}.wait")
    output:
        check_file = DIR_PROC.joinpath("assemblies/verkko/{sample}.{phasing_state}.ok")
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
        assert sorted(processed_partitions) == expected_partitions

        with open(output.check_file, "w") as dump:
            _ = dump.write(f"Processed contigs: {num_contigs}\n")
    # END OF RUN BLOCK


localrules: collect_verkko_output_files
rule collect_verkko_output_files:
    input:
        check_file = DIR_PROC.joinpath("assemblies/verkko/{sample}.{phasing_state}.ok")
    output:
        file_collection = DIR_PROC.joinpath("assemblies/verkko/{sample}.{phasing_state}.output.json")
    conda:
        DIR_ENVS.joinpath("pygraph.yaml")
    params:
        script = find_script("collect_verkko_output"),
        verkko_wd = lambda wildcards, input: pathlib.Path(input.check_file).with_suffix(".wd"),
        base_dir = WORKDIR
    shell:
        "{params.script} --verkko-wd {params.verkko_wd} --base-dir {params.base_dir} "
            "--delete-logs --sample {wildcards.sample} --output {output.file_collection}"


rule compute_verkko_assembly_stats:
    input:
        file_collection = DIR_PROC.joinpath("assemblies/verkko/{sample}.{phasing_state}.output.json")
    output:
        stats = DIR_RES.joinpath(
            "statistics", "assemblies", "{sample}.{phasing_state}.verkko-asm-{asmtype}.statistics.tsv.gz"
        ),
        summary = DIR_RES.joinpath(
            "statistics", "assemblies", "{sample}.{phasing_state}.verkko-asm-{asmtype}.summary.tsv"
        )
    benchmark:
        DIR_RSRC.joinpath("statistics", "assemblies", "{sample}.{phasing_state}.verkko-asm-{asmtype}.stats.rsrc")
    wildcard_constraints:
        asmtype = "(wg|hap1|hap2)"
    conda:
        DIR_ENVS.joinpath("pystats.yaml")
    threads: CPU_HIGH
    resources:
        mem_mb=lambda wildcards, attempt: 24576 * attempt,
        time_hrs=lambda wildcards, attempt: 47 * attempt
    params:
        script=find_script("seqstats"),
        report_seq_lens = " ".join(map(str, [int(1e5), int(5e5), int(1e6), int(1e7), int(5e7), int(1e8)])),
        assembly=lambda wildcards, input: get_verkko_output(
            input.file_collection, f"{wildcards.asmtype}_fasta"
        ),
        acc_res=lambda wildcards, output: register_result(output.stats, output.summary)
    shell:
        "{params.script} --cores {threads} "
        "--summary-length-thresholds {params.report_seq_lens} "
        "--no-homopolymer-runs "
        "--no-canonical-sequence "
        "--temp-records 100 "
        "--str-motif-lengths 2 3 "
        "--output-statistics {output.stats} "
        "--output-summary {output.summary} "
        "--input-files {params.assembly}"


rule compute_verkko_disconn_unassgn_stats:
    input:
        file_collection = DIR_PROC.joinpath("assemblies/verkko/{sample}.{phasing_state}.output.json")
    output:
        stats = DIR_RES.joinpath(
            "statistics", "assemblies", "{sample}.{phasing_state}.verkko-{asmtype}.statistics.tsv.gz"
        ),
        summary = DIR_RES.joinpath(
            "statistics", "assemblies", "{sample}.{phasing_state}.verkko-{asmtype}.summary.tsv"
        )
    benchmark:
        DIR_RSRC.joinpath("statistics", "assemblies", "{sample}.{phasing_state}.verkko-{asmtype}.stats.rsrc")
    wildcard_constraints:
        asmtype = ("disconnected|unassigned")
    conda:
        DIR_ENVS.joinpath("pystats.yaml")
    threads: CPU_MEDIUM
    resources:
        mem_mb=lambda wildcards, attempt: 16384 * attempt,
        time_hrs=lambda wildcards, attempt: 23 * attempt
    params:
        script=find_script("seqstats"),
        report_seq_lens = " ".join(map(str, [int(5e4), int(1e5), int(5e5), int(1e6)])),
        disconn=lambda wildcards, input: get_verkko_output(
            input.file_collection, f"{wildcards.asmtype}_fasta"
        ),
        acc_res=lambda wildcards, output: register_result(output.stats, output.summary)
    shell:
        "{params.script} --cores {threads} "
        "--summary-length-thresholds {params.report_seq_lens} "
        "--no-homopolymer-runs "
        "--no-canonical-sequence "
        "--temp-records 200 "
        "--str-motif-lengths 2 3 "
        "--output-statistics {output.stats} "
        "--output-summary {output.summary} "
        "--input-files {params.disconn}"


localrules: merge_verkko_trio_output
rule merge_verkko_trio_output:
    input:
        vrk = DIR_PROC.joinpath("assemblies/verkko/{sample}.ps-trio")
    output:
        paths = DIR_RES.joinpath(
            "assemblies/verkko/{sample}/{sample}.ps-trio.rukki-paths.tsv"
        ),
        table = DIR_RES.joinpath(
            "assemblies/verkko/{sample}/{sample}.ps-trio.graph-info.tsv"
        ),
    conda:
        DIR_ENVS.joinpath("pygraph.yaml")
    params:
        script = find_script("merge_verkko_infos"),
        graph = lambda wc, input: pathlib.Path(input.vrk).joinpath("assembly.homopolymer-compressed.noseq.gfa"),
        fasta = lambda wc, input: pathlib.Path(input.vrk).joinpath("assembly.fasta"),
        layout = lambda wc, input: pathlib.Path(input.vrk).joinpath("6-layoutContigs/unitig-popped.layout.scfmap"),
        hifi_cov = lambda wc, input: pathlib.Path(input.vrk).joinpath("assembly.hifi-coverage.csv"),
        ont_cov = lambda wc, input: pathlib.Path(input.vrk).joinpath("assembly.ont-coverage.csv"),
        rukki_colors = lambda wc, input: pathlib.Path(input.vrk).joinpath("6-rukki/unitig-popped-unitig-normal-connected-tip.colors.csv"),
        rukki_paths = lambda wc, input: pathlib.Path(input.vrk).joinpath("6-rukki/unitig-popped-unitig-normal-connected-tip.paths.tsv"),
        node_assign = lambda wc, input: pathlib.Path(input.vrk).joinpath("6-rukki/out_final_ann.csv"),
        acc_res = lambda wc, output: register_result(output)
    shell:
        "{params.script} --graph {params.graph} "
        "--fasta {params.fasta} "
        "--scf-layout {params.layout} "
        "--hifi-node-cov {params.hifi_cov} "
        "--ont-node-cov {params.ont_cov} "
        "--rukki-colors {params.rukki_colors} "
        "--rukki-paths {params.rukki_paths} "
        "--rukki-final {params.node_assign} "
        "--out-table {output.table} "
        "--out-path-ids {output.paths}"


localrules: merge_verkko_unphased_output
rule merge_verkko_unphased_output:
    input:
        vrk = DIR_PROC.joinpath("assemblies/verkko/{sample}.ps-none")
    output:
        table = DIR_RES.joinpath(
            "assemblies/verkko/{sample}/{sample}.ps-none.graph-info.tsv"
        ),
    conda:
        DIR_ENVS.joinpath("pygraph.yaml")
    params:
        script = find_script("merge_verkko_infos"),
        graph = lambda wc, input: pathlib.Path(input.vrk).joinpath("assembly.homopolymer-compressed.noseq.gfa"),
        fasta = lambda wc, input: pathlib.Path(input.vrk).joinpath("assembly.fasta"),
        layout = lambda wc, input: pathlib.Path(input.vrk).joinpath("6-layoutContigs/unitig-popped.layout.scfmap"),
        hifi_cov = lambda wc, input: pathlib.Path(input.vrk).joinpath("assembly.hifi-coverage.csv"),
        ont_cov = lambda wc, input: pathlib.Path(input.vrk).joinpath("assembly.ont-coverage.csv"),
        acc_res = lambda wc, output: register_result(output)
    shell:
        "{params.script} --graph {params.graph} "
        "--fasta {params.fasta} "
        "--scf-layout {params.layout} "
        "--hifi-node-cov {params.hifi_cov} "
        "--ont-node-cov {params.ont_cov} "
        "--out-table {output.table} "


rule get_verkko_trio_output_files:
    input:
        file_collect = expand(
            DIR_PROC.joinpath(
                "assemblies/verkko/{sample}.ps-trio.output.json"
            ),
            sample=TRIO_SAMPLES
        ),


rule get_verkko_unphased_output_files:
    input:
        file_collect = expand(
            DIR_PROC.joinpath(
                "assemblies/verkko/{sample}.ps-none.output.json"
            ),
            sample=UNPHASED_SAMPLES
        ),


rule get_verkko_sseq_phased_output_files:
    input:
        file_collect = expand(
            DIR_PROC.joinpath(
                "assemblies/verkko/{sample}.ps-sseq.output.json"
            ),
            sample=SSEQ_SAMPLES
        ),


rule get_verkko_unphased_output_stats:
    input:
        summary = expand(
            DIR_RES.joinpath(
            "statistics", "assemblies", "{sample}.ps-none.verkko-{vrk_out}.summary.tsv"
            ),
            sample=UNPHASED_SAMPLES,
            vrk_out=["asm-wg", "disconnected"]
        ),


rule get_verkko_sseq_phased_output_stats:
    input:
        summary = expand(
            DIR_RES.joinpath(
            "statistics", "assemblies", "{sample}.ps-sseq.verkko-{vrk_out}.summary.tsv"
            ),
            sample=SSEQ_SAMPLES,
            vrk_out=["asm-wg", "asm-hap1", "asm-hap2", "disconnected", "unassigned"]
        ),
