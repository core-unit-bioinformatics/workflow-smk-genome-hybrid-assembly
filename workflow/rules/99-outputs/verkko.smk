
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
        rules.run_verkko_trio_samples.input.assemblies
    )
    VERKKO_OUTPUT.extend(
        rules.get_verkko_trio_graph_info.input.graph_tables
    )

if RUN_VERKKO_UNPHASED_SAMPLES:
    VERKKO_OUTPUT.extend(
        rules.run_verkko_unphased_samples.input.assemblies
    )
    VERKKO_OUTPUT.extend(
        rules.get_verkko_unphased_graph_info.input.graph_tables
    )
