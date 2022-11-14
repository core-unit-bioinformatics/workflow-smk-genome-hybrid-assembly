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
        wd = directory(
            DIR_PROC.joinpath("assemblies/verkko/{sample}.ps-trio")
        ),
        done = DIR_PROC.joinpath("assemblies/verkko/{sample}.ps-trio.ok")
    log:
        DIR_LOG.joinpath("assemblies/verkko/{sample}.ps-trio.log")
    benchmark:
        DIR_RSRC.joinpath("assemblies/verkko/{sample}.ps-trio.rsrc")
    conda:
        "../../envs/verkko.yaml"
    wildcard_constraints:
        sample = CONSTRAINT_TRIO_SAMPLES
    params:
        dryrun = "--dry-run" if VERKKO_DRY_RUN else "",
        check = lambda wildcards, output: "" if VERKKO_DRY_RUN else f" && touch {output.done}",
        acc_in=lambda wildcards, input: register_input(input.hifi, input.nano),
    shell:
        "/usr/bin/time -v "
        "verkko --lsf "
            "--python `which python` "
            "--mbg `which MBG` "
            "--graphaligner `which GraphAligner` "
            "--hifi {input.hifi} "
            "--nano {input.nano} "
            "-d {output.wd} "
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
        wd = directory(
            DIR_PROC.joinpath("assemblies/verkko/{sample}.ps-none")
        ),
        done = DIR_PROC.joinpath("assemblies/verkko/{sample}.ps-none.ok")
    log:
        DIR_LOG.joinpath("assemblies/verkko/{sample}.ps-none.log")
    benchmark:
        DIR_RSRC.joinpath("assemblies/verkko/{sample}.ps-none.rsrc")
    conda:
        "../../envs/verkko.yaml"
    wildcard_constraints:
        sample = CONSTRAINT_UNPHASED_SAMPLES
    params:
        dryrun = "--dry-run" if VERKKO_DRY_RUN else "",
        check = lambda wildcards, output: "" if VERKKO_DRY_RUN else f" && touch {output.done}",
        acc_in=lambda wildcards, input: register_input(input.hifi, input.nano),
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


rule run_verkko_trio_samples:
    input:
        assemblies = expand(
            DIR_PROC.joinpath("assemblies/verkko/{sample}.ps-trio.ok"),
            sample=TRIO_SAMPLES
        ),


rule run_verkko_unphased_samples:
    input:
        assemblies = expand(
            DIR_PROC.joinpath("assemblies/verkko/{sample}.ps-none.ok"),
            sample=UNPHASED_SAMPLES
        ),

