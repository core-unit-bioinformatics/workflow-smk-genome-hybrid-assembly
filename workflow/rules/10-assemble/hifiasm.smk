
rule hifiasm_unphased_samples:
    """
    2023-07-05
    Removed options for debugging
    purposes:
     --write-paf --write-ec
    """
    input:
        hifi = lambda wildcards: MAP_SAMPLE_TO_INPUT_FILES[wildcards.sample]["hifi"],
        nano = lambda wildcards: MAP_SAMPLE_TO_INPUT_FILES[wildcards.sample]["ont"],
    output:
        done = DIR_PROC.joinpath("assemblies/hifiasm/{sample}.ps-none.ok")
    log:
        DIR_LOG.joinpath("assemblies/hifiasm/{sample}.ps-none.log")
    benchmark:
        DIR_RSRC.joinpath("assemblies/hifiasm/{sample}.ps-none.rsrc")
    conda:
        DIR_ENVS.joinpath("hifiasm.yaml")
    wildcard_constraints:
        sample = CONSTRAINT_UNPHASED_SAMPLES
    threads: CPU_MAX
    resources:
        mem_mb=lambda wildcards, attempt: (864 + 576 * attempt) * 1024,
        time_hrs=lambda wildcards, attempt: 167,
    params:
        acc_in=lambda wildcards, input: register_input(input.hifi, input.nano),
        prefix=lambda wildcards, output: pathlib.Path(output.done).with_suffix(".wd").joinpath(wildcards.sample),
        outdir=lambda wildcards, output: pathlib.Path(output.done).with_suffix(".wd")
    shell:
        "mkdir -p {params.outdir}"
            " && "
        "hifiasm -o {params.prefix} -t {threads} "
            "--ul {input.nano} --write-ec "
            "{input.hifi} &> {log}"
            " && "
        "touch {output}"



rule run_hifiasm_unphased_samples:
    input:
        assemblies = expand(
            DIR_PROC.joinpath("assemblies/hifiasm/{sample}.ps-none.ok"),
            sample=UNPHASED_SAMPLES
        ),
