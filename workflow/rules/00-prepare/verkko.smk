localrules: load_verkko_testdata
rule load_verkko_testdata:
    output:
        hifi = DIR_LOCAL_REF.joinpath("hifi.fastq.gz"),
        nano = DIR_LOCAL_REF.joinpath("ont.fastq.gz"),
    shell:
        "curl -L https://obj.umiacs.umd.edu/sergek/shared/ecoli_hifi_subset24x.fastq.gz -o {output.hifi}"
            " && "
        "curl -L https://obj.umiacs.umd.edu/sergek/shared/ecoli_ont_subset50x.fastq.gz -o {output.nano}"


rule run_verkko_test_local:
    input:
        hifi = rules.load_verkko_testdata.output.hifi,
        nano = rules.load_verkko_testdata.output.nano,
    output:
        wd = directory(
            DIR_PROC.joinpath("testdata/verkko/local/assembly")
        ),
        done = DIR_PROC.joinpath("testdata/verkko/local/assembly.ok")
    log:
        DIR_LOG.joinpath("testdata/verkko/local/assembly.log")
    benchmark:
        DIR_RSRC.joinpath("testdata/verkko/local/assembly.log")
    conda:
        DIR_ENVS.joinpath("verkko.yaml")
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: 16384 * attempt,
        mem_gb = lambda wildcards, attempt: 16 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt,
    params:
        dryrun = "--dry-run" if VERKKO_DRY_RUN else "",
        check = lambda wildcards, output: "" if VERKKO_DRY_RUN else f" && touch {output.done}"
    shell:
        "/usr/bin/time -v "
        "verkko --local "
            "--local-memory {resources.mem_gb} "
            "--local-cpus {threads} "
            "--threads {threads} "
            "--python `which python` "
            "--mbg `which MBG` "
            "--graphaligner `which GraphAligner` "
            "--hifi {input.hifi} "
            "--nano {input.nano} "
            "-d {output.wd} "
            "--snakeopts \"--restart-times 1 {params.dryrun}\" "
            " &> {log} {params.check}"


localrules: run_verkko_test_cluster
rule run_verkko_test_cluster:
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
        hifi = rules.load_verkko_testdata.output.hifi,
        nano = rules.load_verkko_testdata.output.nano,
        profile = ancient(config["verkko_smk_profile"]),
    output:
        wd = directory(
            DIR_PROC.joinpath("testdata/verkko/cluster/assembly")
        ),
        done = DIR_PROC.joinpath("testdata/verkko/cluster/assembly.ok")
    log:
        DIR_LOG.joinpath("testdata/verkko/cluster/assembly.log")
    benchmark:
        DIR_RSRC.joinpath("testdata/verkko/cluster/assembly.log")
    conda:
        DIR_ENVS.joinpath("verkko.yaml")
    params:
        dryrun = "--dry-run" if VERKKO_DRY_RUN else "",
        check = lambda wildcards, output: "" if VERKKO_DRY_RUN else f" && touch {output.done}"
    shell:
        "/usr/bin/time -v "
        "verkko --lsf "
            "--python `which python` "
            "--mbg `which MBG` "
            "--graphaligner `which GraphAligner` "
            "--hifi {input.hifi} "
            "--nano {input.nano} "
            "-d {output.wd} "
            "--snakeopts \"--profile $PWD/{input.profile} {params.dryrun}\" "
            "&> {log} {params.check}"


localrules: run_verkko_tests
rule run_verkko_tests:
    input:
        assm_dirs = [
            DIR_PROC.joinpath("testdata/verkko/local/assembly"),
            DIR_PROC.joinpath("testdata/verkko/cluster/assembly")
        ]
