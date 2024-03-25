
READ_STATS = []

if RUN_COMPUTE_ALL_READ_STATS:
    READ_STATS.extend(
        rule.compute_all_read_stats.input.stats
    )
