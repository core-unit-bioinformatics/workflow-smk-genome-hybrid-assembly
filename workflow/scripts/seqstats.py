#!/usr/bin/env python3
import argparse as argp
import collections as col
import functools
import hashlib as hl
import itertools
import multiprocessing as mp
import pathlib as pl
import re
import sys
import time

import dnaio
import numpy as np
import pandas as pd

__prog__ = "seqstats.py"
__version__ = "0.0.1"

# these global consts were introduced
# to avoid using functools.partial
# objects to initialize the sequence
# processors (extra function call = slow)
CONST_ALPHABET = None
CONST_COMPLEMENT_TABLE = None
CONST_LUT_STR_MOTIFS = None
CONST_SEQ_PROCESSORS = None


def parse_command_line():

    parser = argp.ArgumentParser(prog=f"{__prog__} v{__version__}")

    parser.add_argument(
        "--version",
        "-v",
        action="version",
        version=__version__,
        help="Show version and exit."
    )

    parser.add_argument(
        "--n-cpus",
        "--cores",
        "-n",
        type=int,
        default=1,
        dest="cores",
        help="Number of CPU cores to use. Default: 1"
    )

    parser.add_argument(
        "--input-files",
        "--input",
        "-i",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="input_files",
        nargs="+",
        required=True,
        help="Path to input file(s). Also accepts file of filenames (*.fofn)"
    )

    parser.add_argument(
        "--alphabet",
        "-a",
        type=str,
        nargs="+",
        default=list("ACGTN"),
        help="Specify characters of alphabet. Default: A C G T N"
    )

    complement_group = parser.add_mutually_exclusive_group()
    complement_group.add_argument(
        "--complement",
        "-c",
        type=str,
        nargs="+",
        default=list("TGCAN"),
        dest="complement",
        help="Specify complementary characters of alphabet: Default: T G C A N"
    )
    complement_group.add_argument(
        "--no-canonical-sequence",
        "-ncan",
        action="store_true",
        default=False,
        dest="no_canonical_sequence",
        help="If set, do not use the canonical sequence for computations. Default: False"
    )

    parser.add_argument(
        "--no-sequence-composition",
        "-nsc",
        action="store_true",
        default=False,
        dest="no_sequence_composition",
        help="If set, do not compute sequence composition. Default: False"
    )

    parser.add_argument(
        "--no-homopolymer-runs",
        "-nhpr",
        action="store_true",
        default=False,
        dest="no_homopolymer_runs",
        help="If set, do not count homopolymer runs in sequence. Default: False"
    )

    str_group = parser.add_mutually_exclusive_group()
    str_group.add_argument(
        "--no-short-tandem-repeats",
        "-nstr",
        action="store_true",
        default=False,
        dest="no_short_tandem_repeats",
        help="If set, do not count short tandem repeats (motif length < 7) in sequence. Default: False"
    )

    str_group.add_argument(
        "--str-motif-lengths",
        "-strl",
        type=int,
        nargs="+",
        default=[2, 3],
        dest="str_motif_lengths",
        help="Specify motif length for STR counting. Default: 2 3 4"
    )

    parser.add_argument(
        "--summary-length-thresholds",
        "-slt",
        type=int,
        nargs="+",
        dest="length_thresholds",
        default=[20000, 50000, 100000]
    )

    parser.add_argument(
        "--coverage-ref-size",
        "-crs",
        type=int,
        default=3000000000,
        dest="coverage_ref_size"
    )

    args = parser.parse_args()
    return args


def seq_norm_noncanonical(sequence):
    return sequence.upper()


def seq_norm_canonical(sequence):

    sequence = sequence.upper()
    rc_seq = revcomp(sequence)
    if sequence <= rc_seq:
        return sequence
    return rc_seq


def revcomp(sequence):
    global CONST_COMPLEMENT_TABLE
    return sequence.translate(CONST_COMPLEMENT_TABLE)[::-1]


def set_alphabet_complement_table(complement):

    if len(CONST_ALPHABET) != len(complement):
        raise ValueError(
            "Alphabet and complement must be of same length:\n"
            f"{CONST_ALPHABET}\n"
            f"{complement}\n"
        )
    complement_table = str.maketrans(
        dict((a, c) for a, c in zip(CONST_ALPHABET, complement))
    )

    global CONST_COMPLEMENT_TABLE
    CONST_COMPLEMENT_TABLE = complement_table

    return


def sequence_composition(sequence, seq_info):
    char_counts = col.Counter(sequence)
    for k, v in char_counts.items():
        seq_info[f"char{k}_cov"] = v
    return


def count_homopolymer_runs(sequence, seq_info):

    distinct = 0
    total_cov = 0
    run_lengths = col.defaultdict(list)
    for nuc in CONST_ALPHABET:
        pattern = f'{nuc}{{2,}}'
        for mobj in re.finditer(pattern, sequence):
            s, e = mobj.span()
            cov = e - s
            assert cov > 1
            seq_info[f"hpr{nuc}_cov"] += cov
            run_lengths[nuc].append(cov)
            total_cov += cov
            distinct += 1
    seq_info["hpr_cov"] = cov
    seq_info["hpr_length"] = len(sequence) - total_cov + distinct
    for nuc, nuc_hp_rl in run_lengths.items():
        seq_info[f"hpr{nuc}_rl_max"] = max(nuc_hp_rl)
        seq_info[f"hpr{nuc}_rl_median"] = sorted(nuc_hp_rl)[len(nuc_hp_rl)//2]
    return


def count_short_tandem_repeats(motifs, sequence, seq_info):

    total_cov = 0
    repeat_nums = col.defaultdict(list)
    motif_length = len(motifs[0])
    for motif in motifs:
        pattern = f'({motif}){{2,}}'
        for mobj in re.finditer(pattern, sequence):
            s, e = mobj.span()
            cov = e - s
            assert cov > motif_length
            seq_info[f'rep{motif_length}_{motif}_cov'] += cov
            repeat_nums[motif].append(cov/motif_length)
            total_cov += cov
    seq_info[f'rep{motif_length}_total_cov'] = cov
    for motif, motif_repeat_nums in repeat_nums.items():
        seq_info[f'rep{motif_length}_{motif}_rn_max'] = max(motif_repeat_nums)
        seq_info[f'rep{motif_length}_{motif}_rn_median'] = sorted(motif_repeat_nums)[len(motif_repeat_nums)//2]
    return


def get_repeat_motifs(motif_length):

    motifs = []
    for repeat_motif in itertools.product(CONST_ALPHABET, repeat=motif_length):
        if len(set(repeat_motif)) == 1:
            continue
        motifs.append(''.join(repeat_motif))
    return motifs


def process_sequence(seq_record):

    start = time.perf_counter()

    make_canonical, _ = CONST_SEQ_PROCESSORS[0]
    if make_canonical:
        canon_seq = seq_norm_canonical(seq_record.sequence)
    else:
        canon_seq = seq_norm_noncanonical(seq_record.sequence)

    seq_hash = hl.sha256(canon_seq.encode("utf-8")).hexdigest()

    seq_stats = col.Counter()
    seq_stats["seq_length"] = len(canon_seq)
    if len(CONST_SEQ_PROCESSORS) > 1:
        for (processor, proc_spec) in CONST_SEQ_PROCESSORS[1:]:
            if processor == "sequence_composition":
                sequence_composition(canon_seq, seq_stats)
            elif processor == "homopolymer_runs":
                count_homopolymer_runs(canon_seq, seq_stats)
            elif processor == "short_tandem_repeats":
                motifs = CONST_LUT_STR_MOTIFS[proc_spec]
                count_short_tandem_repeats(motifs, canon_seq, seq_stats)
            else:
                raise ValueError(f"Unknown sequence processor: {processor}")

    end = time.perf_counter()
    proc_time = end-start

    return seq_record.name, seq_hash, seq_stats, proc_time


def setup_sequence_processors(args):

    processors = []
    if args.no_canonical_sequence:
        seqnorm = ("noncanonical", None)
    else:
        seqnorm = ("canonical", None)
    processors.append(seqnorm)

    if not args.no_sequence_composition:
        processors.append(("sequence_composition", None))

    if not args.no_homopolymer_runs:
        processors.append(("homopolymer_runs", None))

    if not args.no_short_tandem_repeats:
        global CONST_LUT_STR_MOTIFS
        CONST_LUT_STR_MOTIFS = dict()
        for motif_length in args.str_motif_lengths:
            motifs = get_repeat_motifs(motif_length)
            CONST_LUT_STR_MOTIFS[motif_length] = motifs
            processors.append(("short_tandem_repeats", motif_length))

    global CONST_SEQ_PROCESSORS
    CONST_SEQ_PROCESSORS = processors

    return


def compute_length_statistics(seq_source, seq_lengths, threshold, ref_size):

    tmp_lengths = np.sort(seq_lengths[seq_lengths > threshold])[::-1]
    length_stats = []
    if tmp_lengths.size > 0:
        total_length = tmp_lengths.sum()
        # NB: sorted in descending order
        cum_length = tmp_lengths.cumsum()

        length_stats.append(
            (seq_source, f"total_length_grt_{threshold}", total_length)
        )

        length_stats.append(
            (seq_source, f"total_length_grt_{threshold}", tmp_lengths.size)
        )

        length_stats.append(
            (seq_source, f"cov_xfold_at_{ref_size}",
             round(total_length / ref_size, 1))
        )

        length_stats.append(
            (seq_source, f"length_N50_grt_{threshold}",
             tmp_lengths[cum_length>(total_length//2)].max())
        )

        length_stats.append(
            (seq_source, f"length_auN_grt_{threshold}",
             int(round(sum(l*l/total_length for l in tmp_lengths), 0)))
        )

    return length_stats


def prepare_summary(args, stats, proc_timings):

    proc_timings = np.sort(np.array(proc_timings, dtype=float))
    proc_median = proc_timings[proc_timings.size//2]
    proc_mean = proc_timings.mean()
    proc_top = proc_timings[int(proc_timings.size * 0.99)]

    ref_size = args.coverage_ref_size
    length_thresholds = args.length_thresholds
    if min(length_thresholds) != 0:
        length_thresholds = [0] + length_thresholds

    summary_stats = [
        ("all", "sec_per_read_median", proc_median),
        ("all", "sec_per_read_mean", proc_mean),
        ("all", "sec_per_read_99pct", proc_top)
    ]
    for t in length_thresholds:
        t_len_stats = compute_length_statistics(
            "all", stats["seq_length"].values, t, ref_size
        )
        summary_stats.extend(t_len_stats)

    seq_sources = stats.index.unique(level="seq_source")
    if seq_sources.size > 1:

        for seq_source in seq_sources:
            sub = stats.xs(seq_source, level="seq_source")
            for t in length_thresholds:
                t_len_stats = compute_length_statistics(
                    seq_source, sub["seq_length"].values, t, ref_size
                )
                summary_stats.extend(t_len_stats)

    summary_df = pd.DataFrame.from_records(
        summary_stats,
        columns=["source", "statistic", "value"]
    )
    return summary_df


def main():

    args = parse_command_line()

    global CONST_ALPHABET
    CONST_ALPHABET = args.alphabet

    _ = set_alphabet_complement_table(args.complement)

    all_inputs = args.input_files
    if isinstance(all_inputs, list) and len(all_inputs) == 1:
        if all_inputs[0].suffix in ["fofn", "fof", "fofp"]:
            all_inputs = []
            with open(all_inputs[0], "r") as fofn:
                for line in fofn:
                    this_file = pl.Path(line.strip()).resolve(strict=True)
                    all_inputs.append(this_file)

    setup_sequence_processors(args)

    stats = []
    index_records = []
    proc_timings = []
    with mp.Pool(args.cores) as pool:
        for input_file in sorted(all_inputs):
            file_name = input_file.name
            with dnaio.open(str(input_file)) as seq_file:
                stats_iter = pool.imap_unordered(process_sequence, seq_file)
                for seq_name, seq_hash, seq_stats, proc_time in stats_iter:
                    index_records.append((seq_name, file_name, seq_hash))
                    stats.append(seq_stats)
                    proc_timings.append(proc_time)

    mindex = pd.MultiIndex.from_tuples(index_records, names=["seq_name", "seq_source", "seq_hash"])
    stats = pd.DataFrame.from_records(stats, index=mindex)
    stats.fillna(0, inplace=True)

    summary = prepare_summary(args, stats, proc_timings)
    print(summary)

    return 0



if __name__ == "__main__":
    main()
