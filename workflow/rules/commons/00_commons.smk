# some constants built from the Snakemake
# command line arguments
include: "01_constants.smk"
# Module containing generic
# Python utility functions
include: "02_pyutils.smk"
# Module containing generic
# Snakemake utility rules, e.g., dump_config
include: "03_smkutils.smk"
# Module containing
# reference container location information
include: "05_refcon.smk"
# Module performing state/env-altering
# operations before Snakemake starts its
# actual work
include: "09_staging.smk"
