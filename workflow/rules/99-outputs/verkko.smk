
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
        rules.get_verkko_trio_output_files.input.file_collect
    )

if RUN_VERKKO_UNPHASED_SAMPLES:
    VERKKO_OUTPUT.extend(
        rules.get_verkko_unphased_output_files.input.file_collect
    )
    VERKKO_OUTPUT.extend(
        rules.get_verkko_unphased_output_stats.input.summary
    )

if RUN_VERKKO_SSEQ_SAMPLES:
    VERKKO_OUTPUT.extend(
        rules.get_verkko_sseq_phased_output_files.input.file_collect
    )
    VERKKO_OUTPUT.extend(
        rules.get_verkko_sseq_phased_output_stats.input.summary
    )
