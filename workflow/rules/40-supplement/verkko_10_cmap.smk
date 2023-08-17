
rule homopolymer_compress_verkko_whole_genome:
    input:
        fasta = expand(
            rules.filter_verkko_dup_sequences.output.asm_unit,
            asm_unit=["hap1", "hap2", "unassigned", "disconnected", "ebv", "mito", "rdna"],
            allow_missing=True
        ),
        faidx = expand(
            rules.filter_verkko_dup_sequences.output.fai,
            asm_unit=["hap1", "hap2", "unassigned", "disconnected", "ebv", "mito", "rdna"],
            allow_missing=True
        )
    output:
        fasta = DIR_PROC.joinpath(
            "40-supplement", "verkko", "fasta_seq",
            "{sample}.{phasing_state}.fastaseq.hpc.fasta.gz"
        ),
        faidx = DIR_PROC.joinpath(
            "40-supplement", "verkko", "fasta_seq",
            "{sample}.{phasing_state}.fastaseq.hpc.fasta.gz.fai"
        ),
        cmap = DIR_PROC.joinpath(
            "40-supplement", "verkko", "fasta_seq",
            "{sample}.{phasing_state}.cmap.hpc.tsv.gz"
        ),
        tbi = DIR_PROC.joinpath(
            "40-supplement", "verkko", "fasta_seq",
            "{sample}.{phasing_state}.cmap.hpc.tsv.gz.tbi"
        ),
    benchmark:
        DIR_RSRC.joinpath(
            "40-supplement", "verkko", "fasta_seq",
            "{sample}.{phasing_state}.fastaseq.hpc.rsrc"
        ),
    log:
        DIR_LOG.joinpath(
            "40-supplement", "verkko", "fasta_seq",
            "{sample}.{phasing_state}.fastaseq.hpc.log"
        ),
    conda: DIR_ENVS.joinpath("pyseq.yaml")
    threads: CPU_LOW
    resources:
        mem_mb = lambda wildcards, attempt: 16384 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script = find_script("seq_hpc"),
        plain_tsv = lambda wildcards, output: pathlib.Path(output.cmap).with_suffix("")
    shell:
        "zcat {input.fasta} | {params.script} --cmap-table {params.plain_tsv} --report "
        "| bgzip -c > {output.fasta}"
            " && "
        "samtools faidx {output.fasta}"
            " && "
        "bgzip --threads {threads} {params.plain_tsv}"
            " && "
        "tabix --zero-based --sequence 1 --begin 2 --end 3 --comment \"#\" {output.cmap}"


rule minimap_align_verkko_graphseq_to_fastaseq:
    input:
        fastaseq = rules.homopolymer_compress_verkko_whole_genome.output.fasta,
        gfaseq = get_verkko_gfaseq_hpc_fasta
    output:
        paf = DIR_PROC.joinpath(
            "40-supplement", "verkko", "gfa_to_fasta_align",
            "{sample}.{phasing_state}.paf.gz"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: CPU_MEDIUM
    resources:
        time_hrs = lambda wildcards, attempt: 1 * attempt,
        mem_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
    shell:
        "minimap2 -c -x asm5 -N 1 -p 0.95 --cs --eqx -t {threads} {input.fastaseq} {input.gfaseq}"
            " | "
        "pigz -p {threads} > {output.paf}"


rule normalize_minimap_gfa_to_fasta_align_paf:
    input:
        paf = rules.minimap_align_verkko_graphseq_to_fastaseq.output.paf
    output:
        tsv = DIR_PROC.joinpath(
            "40-supplement", "verkko", "gfa_to_fasta_align",
            "{sample}.{phasing_state}.norm-paf.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    params:
        script=find_script("normalize_paf")
    shell:
        "{params.script} --input {input.paf} --output {output.tsv}"


rule run_verkko_supplement_cmap:
    input:
        cmap = expand(
            rules.homopolymer_compress_verkko_whole_genome.output.cmap,
            sample=SSEQ_SAMPLES,
            phasing_state=["ps-sseq"]
        ),
        #tsv = expand(
        #    rules.normalize_minimap_gfa_to_fasta_align_paf.output.tsv,
        #    sample=SSEQ_SAMPLES,
        #    phasing_state=["ps-sseq"]
        #)

