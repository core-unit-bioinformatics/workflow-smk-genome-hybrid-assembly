# Documentation for Snakemake workflow "genome hybrid assembly"

**Disclaimer** the workflow is in prototype state and configuration may change at any time

This workflow creates state-of-the-art genome hybrid assemblies
for diploid vertebrate species. This version of the workflow was
developed for the following scenario:

- input species: human
  - successfully tested as well: muntjac
- assembler: Verkko v1.4.1
  - hifiasm v0.19.x is not yet fully integrated
- inputs:
  - long accurate reads: PacBio HiFi (Sequel-II/Revio)
    - required coverage: at least ~40X, ideally ~60X
  - long connecting reads: Oxford Nanopore ultralong (R9)
    - required coverage: ~30X ultralong (>100 kbp) reads
  - optional input for phasing:
    - trio: kmer databases created with meryl
    - HiC: HiC short reads
    - graph/node coloring for Verkko's Rukki (GFA file)
- outputs:
  - main: whole-genome assembly, potentially phased
  - main: basic (length) statistics about assembly and long reads
  - optional: a coordinate map between the homopolymer-compressed assembly graph and the linearized plain FASTA files

The sample sheet must be a TAB-separated text file (`.tsv` file extension) with at least the columns
`sample`, `hifi` and `ont`, where both `hifi` and `ont` columns can hold an arbitrary number of
input file paths (comma-seperated, i.e., `file_path1,file_path2,file_path3`) representing the respective
read dataset for that sample. Common file extensions are recognised (e.g., `fastq.gz`, `.fq.gz` and so on).
The Verkko assembler can optionally be configured for using three different phasing signals;
add the column `target` to the sample sheet plus the following fields:
1. trio-based: set value `trio` in column `target` and add columns `hap1` and `hap2` pointing to meryl k-mer databases
of the sample parents (conventionally, `hap1` should be the father and `hap2` the mother)
2. Hi-C: set value `hic` in column `target` and add fields `hic1` and `hic2` for the Hi-C reads of mate 1 and 2, respectively
3. Strand-seq: set value `sseq` in column `target` and add field `phasing_paths` pointing to a `.gaf` format file produced
by the [Grapahasing pipeline](https://github.com/marschall-lab/strand-seq-graph-phasing)

Since Verkko itself is implemented as a Snakemake workflow, you can execute a dry run to check if all
input requirements are met by setting the option `verkko_dry_run` to `true`, see this example configuration:

[Example parameterization](../config/parameter.yaml)

## User documentation for workflow template

All standard workflows of the CUBI implement the same user
interface (or at least aim for a highly similar interface).
Hence, before [executing the workflow](concepts/running.md),
we strongly recommend reading through the documentation
that explains how we help you to keep track of your analysis
results; we refer to this concept as
[**"file accounting"**](concepts/accounting.md). This feature
of standard CUBI workflows enables the pipeline to auto-
matically create a so-called [**"manifest"** file](concepts/accounting.md)
for your analysis run.

In case of questions, please open a GitHub issue in the repository
of the workflow you are trying to execute.

## Developer documentation

Besides reading the user documentation, CUBI developers find more
information regarding standadized workflow development in the
[developer notes](concepts/developing.md). Please keep in mind
to always cross-link that information with the guidelines
published in the
[CUBI knowledge base](https://github.com/core-unit-bioinformatics/knowledge-base/wiki/).

Please raise any issues with these guidelines "close to the code",
i.e., either open an issue in the
[knowledge base repo](https://github.com/core-unit-bioinformatics/knowledge-base)
or in the affected repo for more specific cases.
