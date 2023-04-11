"""
Use this module to list all includes
required for your pipeline - do not
add your pipeline-specific modules
to "commons/00_commons.smk"
"""

include: "00-prepare/settings.smk"
include: "00-prepare/sample_table.smk"
include: "00-prepare/verkko.smk"

include: "10-assemble/verkko_pyutils.smk"
include: "10-assemble/verkko.smk"
include: "30-postprocess/verkko.smk"

include: "99-outputs/verkko.smk"
