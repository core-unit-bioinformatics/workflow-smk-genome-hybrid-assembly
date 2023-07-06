
DATA_ROOT = config.get("data_root", "/")
DATA_ROOT = pathlib.Path(DATA_ROOT).resolve(strict=True)

VERKKO_DRY_RUN = config.get("verkko_dry_run", False)

RUN_VERKKO_TEST_LOCAL = config.get("run_verkko_test_local", True)
RUN_VERKKO_TEST_CLUSTER = config.get("run_verkko_test_cluster", True)

RUN_VERKKO_TRIO_SAMPLES = config.get("run_verkko_trio_samples", True)
RUN_VERKKO_UNPHASED_SAMPLES = config.get("run_verkko_unphased_samples", True)
RUN_VERKKO_SSEQ_SAMPLES = config.get("run_verkko_sseq_samples", True)
