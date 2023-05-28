
rule compute_input_stats:
    input:
        reads=lambda wildcards: MAP_SAMPLE_TO_INPUT_FILES[wildcards.sample][wildcards.read_type]
    output:
        stats = DIR_RES.joinpath(
            "statistics", "input", "{sample}_{read_type}.statistics.tsv.gz"
        ),
        summary = DIR_RES.joinpath(
            "statistics", "input", "{sample}_{read_type}.summary.tsv"
        )
    benchmark:
        DIR_RSRC.joinpath("statistics", "input", "{sample}_{read_type}.stats.rsrc")
    wildcard_constraints:
        read_type="(hifi|ont)"
    conda:
        DIR_ENVS.joinpath("pystats.yaml")
    threads: CPU_MEDIUM
    resources:
        mem_mb=lambda wildcards, attempt: 32768 * attempt,
        time_hrs=lambda wildcards, attempt: 11*attempt
    params:
        script=find_script("seqstats"),
        report_seq_lens=lambda wildcards: (
            " 10000 15000 20000 50000 " if wildcards.read_type == "hifi"
            else " 25000 50000 100000 500000 1000000 "
        ),
        timings_out=lambda wildcards, output: (
            f" --output-timings {output.stats.replace('summary.tsv', 'proc-timings.tsv.gz')} "
            if "NA" in wildcards.sample else " "
        ),
        acc_res=lambda wildcards, output: register_result(output.stats, output.summary)
    shell:
        "{params.script} --cores {threads} "
        "--summary-length-thresholds {params.report_seq_lens} "
        "--temp-records 100000 "
        "--output-statistics {output.stats} "
        "--output-summary {output.summary} "
        "{params.timings_out}"
        "--input-files {input.reads}"


rule compute_all_input_stats:
    input:
        stats = expand(DIR_RES.joinpath(
                "statistics", "input",
                "{sample}_{read_type}.summary.tsv"
            ),
            sample=SAMPLES,
            read_type=["hifi", "ont"]
        )
