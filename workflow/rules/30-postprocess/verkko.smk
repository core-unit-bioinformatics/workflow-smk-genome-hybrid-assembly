
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
        "../../envs/pygraph.yaml"
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
        "../../envs/pygraph.yaml"
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


rule get_verkko_trio_graph_info:
    input:
        graph_tables = expand(DIR_RES.joinpath(
            "assemblies/verkko/{sample}/{sample}.ps-trio.graph-info.tsv"),
            sample=TRIO_SAMPLES
        ),


rule get_verkko_unphased_graph_info:
    input:
        graph_tables = expand(DIR_RES.joinpath(
            "assemblies/verkko/{sample}/{sample}.ps-none.graph-info.tsv"),
            sample=UNPHASED_SAMPLES
        ),
