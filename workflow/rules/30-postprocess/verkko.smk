
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
            "statistics", "assemblies", "{sample}.{phasing_state}.verkko-asm.statistics.tsv.gz"
        ),
        summary = DIR_RES.joinpath(
            "statistics", "assemblies", "{sample}.{phasing_state}.verkko-asm.summary.tsv"
        )
    benchmark:
        DIR_RSRC.joinpath("statistics", "assemblies", "{sample}.{phasing_state}.verkko-asm.stats.rsrc")
    conda:
        DIR_ENVS.joinpath("pystats.yaml")
    threads: CPU_MEDIUM
    resources:
        mem_mb=lambda wildcards, attempt: 16384 * attempt,
        time_hrs=lambda wildcards, attempt: 23*attempt
    params:
        script=find_script("seqstats"),
        report_seq_lens = " ".join(map(str, [int(1e5), int(5e5), int(1e6), int(1e7), int(5e7), int(1e8)])),
        assembly=lambda wildcards, input: get_verkko_output(input.file_collection, "wg_fasta"),
        acc_res=lambda wildcards, output: register_result(output.stats, output.summary)
    shell:
        "{params.script} --cores {threads} "
        "--summary-length-thresholds {params.report_seq_lens} "
        "--no-homopolymer-runs "
        "--no-canonical-sequence "
        "--temp-records 200 "
        "--str-motif-lengths 2 "
        "--output-statistics {output.stats} "
        "--output-summary {output.summary} "
        "--input-files {params.assembly}"


rule compute_verkko_unassembled_stats:
    input:
        file_collection = DIR_PROC.joinpath("assemblies/verkko/{sample}.{phasing_state}.output.json")
    output:
        stats = DIR_RES.joinpath(
            "statistics", "assemblies", "{sample}.{phasing_state}.verkko-disconn.statistics.tsv.gz"
        ),
        summary = DIR_RES.joinpath(
            "statistics", "assemblies", "{sample}.{phasing_state}.verkko-disconn.summary.tsv"
        )
    benchmark:
        DIR_RSRC.joinpath("statistics", "assemblies", "{sample}.{phasing_state}.verkko-disconn.stats.rsrc")
    conda:
        DIR_ENVS.joinpath("pystats.yaml")
    threads: CPU_LOW
    resources:
        mem_mb=lambda wildcards, attempt: 8192 * attempt,
        time_hrs=lambda wildcards, attempt: 11 * attempt
    params:
        script=find_script("seqstats"),
        report_seq_lens = " ".join(map(str, [int(5e4), int(1e5), int(5e5), int(1e6)])),
        disconn=lambda wildcards, input: get_verkko_output(input.file_collection, "disconnected"),
        acc_res=lambda wildcards, output: register_result(output.stats, output.summary)
    shell:
        "{params.script} --cores {threads} "
        "--summary-length-thresholds {params.report_seq_lens} "
        "--no-homopolymer-runs "
        "--no-canonical-sequence "
        "--temp-records 200 "
        "--str-motif-lengths 2 "
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


rule get_verkko_unphased_output_stats:
    input:
        summary = expand(
            DIR_RES.joinpath(
            "statistics", "assemblies", "{sample}.ps-none.verkko-{vrk_out}.summary.tsv"
            ),
            sample=UNPHASED_SAMPLES,
            vrk_out=["asm", "disconn"]
        ),
