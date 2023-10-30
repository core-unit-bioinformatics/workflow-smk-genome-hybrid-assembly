
VERKKO_OUTPUT = []

if RUN_VERKKO_TEST_LOCAL and RUN_VERKKO_TEST_CLUSTER:
    VERKKO_OUTPUT.extend(
        rules.run_verkko_tests.input.assm_dirs
    )
elif RUN_VERKKO_TEST_LOCAL:
    VERKKO_OUTPUT.extend(
        rules.run_verkko_test_local.output.wd
    )
elif RUN_VERKKO_TEST_CLUSTER:
    VERKKO_OUTPUT.extend(
        rules.run_verkko_test_cluster.output.wd
    )
else:
    pass

if RUN_VERKKO_TRIO_SAMPLES:
    VERKKO_OUTPUT.extend(
        rules.get_verkko_trio_phased_output_stats.input.summary
    )

if RUN_VERKKO_UNPHASED_SAMPLES:
    VERKKO_OUTPUT.extend(
        rules.get_verkko_unphased_output_stats.input.summary
    )

if RUN_VERKKO_SSEQ_SAMPLES:
    VERKKO_OUTPUT.extend(
        rules.postprocess_verkko_sseq_samples.input.exemplars
    )
    VERKKO_OUTPUT.extend(
        rules.postprocess_verkko_sseq_samples.input.asm_units
    )
    VERKKO_OUTPUT.extend(
        rules.get_verkko_sseq_phased_output_stats.input.summary
    )


if RUN_VERKKO_HIC_SAMPLES:
    VERKKO_OUTPUT.extend(
        rules.postprocess_verkko_hic_samples.input.exemplars
    )
    VERKKO_OUTPUT.extend(
        rules.postprocess_verkko_hic_samples.input.asm_units
    )
    VERKKO_OUTPUT.extend(
        rules.get_verkko_hic_phased_output_stats.input.summary
    )
