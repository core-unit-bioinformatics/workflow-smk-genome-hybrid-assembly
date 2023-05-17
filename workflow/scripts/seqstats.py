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
        default=list("ACGT"),
        help="Specify characters of alphabet. Default: A C G T"
    )

    complement_group = parser.add_mutually_exclusive_group()
    complement_group.add_argument(
        "--complement",
        "-c",
        type=str,
        nargs="+",
        default=list("TGCA"),
        dest="complement",
        help="Specify complementary characters of alphabet: Default: T G C A"
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

    parser.add_argument(
        "--output-statistics",
        "-o",
        type=lambda x: pl.Path(x).resolve(),
        default=None,
        dest="output_statistics",
        help="Path to output TSV file to dump statistics. Default: None"
    )

    parser.add_argument(
        "--output-summary",
        "-s",
        type=lambda x: pl.Path(x).resolve(),
        default=None,
        dest="output_summary",
        help="Path to output summary TSV. Default: None"
    )

    parser.add_argument(
        "--output-timings",
        "-t",
        type=lambda x: pl.Path(x).resolve(),
        default=None,
        dest="output_timings",
        help="Path to output processing timings per read. Default: None"
    )

    args = parser.parse_args()

    no_stats = args.output_statistics is None
    no_summary = args.output_summary is None
    no_timings = args.output_timings is None

    if all([no_stats, no_summary, no_timings]):
        raise RuntimeError("Not a single output file specified - aborting.")

    return args


def seq_norm_noncanonical(sequence):
    return sequence.upper()


def seq_norm_canonical(complement_table, sequence):

    sequence = sequence.upper()
    rc_seq = sequence.translate(complement_table)[::-1]
    if sequence <= rc_seq:
        canonical_seq = sequence
    else:
        canonical_seq = rc_seq
    return canonical_seq


def get_alphabet_complement_table(alphabet, complement):

    if len(alphabet) != len(complement):
        raise ValueError(
            "Alphabet and complement must be of same length:\n"
            f"{alphabet}\n"
            f"{complement}\n"
        )
    complement_table = str.maketrans(
        dict((a, c) for a, c in zip(alphabet, complement))
    )

    return complement_table


def sequence_composition(sequence, seq_stats):
    char_counts = col.Counter(sequence)
    for k, v in char_counts.items():
        seq_stats[f"char{k}_cov"] = v
    return


def count_homopolymer_runs(alphabet, sequence, seq_stats):

    distinct = 0
    total_cov = 0
    run_lengths = col.defaultdict(list)
    for nuc in alphabet:
        pattern = f'{nuc}{{2,}}'
        for mobj in re.finditer(pattern, sequence):
            s, e = mobj.span()
            cov = e - s
            assert cov > 1
            seq_stats[f"hpr{nuc}_cov"] += cov
            run_lengths[nuc].append(cov)
            total_cov += cov
            distinct += 1
    seq_stats["hpr_cov"] = cov
    seq_stats["hpr_length"] = len(sequence) - total_cov + distinct
    for nuc, nuc_hp_rl in run_lengths.items():
        seq_stats[f"hpr{nuc}_rl_max"] = max(nuc_hp_rl)
        seq_stats[f"hpr{nuc}_rl_median"] = sorted(nuc_hp_rl)[len(nuc_hp_rl)//2]
    return


def count_short_tandem_repeats(motifs, sequence, seq_stats):

    total_cov = 0
    repeat_nums = col.defaultdict(list)
    motif_length = len(motifs[0])
    for motif in motifs:
        pattern = f'({motif}){{2,}}'
        for mobj in re.finditer(pattern, sequence):
            s, e = mobj.span()
            cov = e - s
            assert cov > motif_length
            seq_stats[f'rep{motif_length}_{motif}_cov'] += cov
            repeat_nums[motif].append(cov/motif_length)
            total_cov += cov
    seq_stats[f'rep{motif_length}_total_cov'] = total_cov
    for motif, motif_repeat_nums in repeat_nums.items():
        seq_stats[f'rep{motif_length}_{motif}_rn_max'] = max(motif_repeat_nums)
        seq_stats[f'rep{motif_length}_{motif}_rn_median'] = sorted(motif_repeat_nums)[len(motif_repeat_nums)//2]
    return


def get_repeat_motifs(alphabet, motif_length):

    motifs = []
    for repeat_motif in itertools.product(alphabet, repeat=motif_length):
        if len(set(repeat_motif)) == 1:
            continue
        motifs.append(''.join(repeat_motif))
    return motifs


def process_sequence(processors, seq_record):

    proc_times = {
        "seq_name": seq_record.name,
        "seq_length": len(seq_record.sequence)
    }

    total_start = time.perf_counter()

    seq_normalizer, proc_name = processors[0]
    norm_start = time.perf_counter()
    norm_seq = seq_normalizer(seq_record.sequence)
    norm_end = time.perf_counter()
    proc_times[proc_name] = norm_end - norm_start

    seq_hash = hl.sha256(norm_seq.encode("utf-8")).hexdigest()

    seq_stats = col.Counter()
    seq_stats["seq_length"] = len(norm_seq)


    for processor, proc_name in processors[1:]:
        proc_start = time.perf_counter()
        processor(norm_seq, seq_stats)
        proc_end = time.perf_counter()
        proc_times[proc_name] = proc_end - proc_start

    total_end = time.perf_counter()
    proc_times["total"] = total_end - total_start

    return seq_record.name, seq_hash, seq_stats, proc_times


def setup_sequence_processors(args):

    processors = []
    if args.no_canonical_sequence:
        seqnorm = seq_norm_noncanonical
    else:
        complement_table = get_alphabet_complement_table(
            args.alphabet, args.complement
        )
        seqnorm = functools.partial(
            seq_norm_canonical,
            complement_table
        )
    processors.append((seqnorm, "seq_normalize"))

    if not args.no_sequence_composition:
        processors.append((sequence_composition, "seq_composition"))

    if not args.no_homopolymer_runs:
        hpr_init = functools.partial(
            count_homopolymer_runs,
            args.alphabet
        )
        processors.append((hpr_init, "homopolymer_runs"))

    if not args.no_short_tandem_repeats:
        for motif_length in args.str_motif_lengths:
            motifs = get_repeat_motifs(args.alphabet, motif_length)
            str_proc_init = functools.partial(
                count_short_tandem_repeats,
                tuple(motifs)
            )
            proc_name = f"tandemrep_len{motif_length}"
            processors.append((str_proc_init, proc_name))

    return processors


def compute_length_statistics(seq_source, seq_lengths, threshold, ref_size):

    tmp_lengths = np.sort(seq_lengths[seq_lengths > threshold])[::-1]
    length_stats = []
    if tmp_lengths.size > 0:

        readable_threshold = make_number_human_readable(threshold)
        readable_refsize = make_number_human_readable(ref_size)

        total_length = int(tmp_lengths.sum())
        # NB: sorted in descending order
        cum_length = tmp_lengths.cumsum()

        length_stats.append(
            (seq_source, f"total_length_grt_{readable_threshold}", total_length)
        )

        length_stats.append(
            (seq_source, f"total_num_grt_{readable_threshold}", tmp_lengths.size)
        )

        length_stats.append(
            (seq_source, f"cov_xfold_at_{readable_refsize}",
             round(total_length / ref_size, 1))
        )

        length_stats.append(
            (seq_source, f"length_N50_grt_{readable_threshold}",
             int(tmp_lengths[cum_length>(total_length//2)].max()))
        )

        length_stats.append(
            (seq_source, f"length_auN_grt_{readable_threshold}",
             int(round(sum(l*l/total_length for l in tmp_lengths), 0)))
        )

    return length_stats


def make_number_human_readable(number):

    readable_number = None

    if number < 1000:
        readable_number = f"{number}bp"
    elif 1000 <= number < int(1e6):
        if number % 1000 != 0:
            readable_number = f"{number}bp"
        else:
            kilo_num = int(number / 1000)
            readable_number = f"{kilo_num}kbp"
    elif int(1e6) <= number < int(1e9):
        if number % int(1e6) != 0:
            readable_number = f"{number}bp"
        else:
            mega_num = int(number / 1e6)
            readable_number = f"{mega_num}Mbp"
    elif int(1e9) <= number < int(1e12):
        if number % int(1e9) != 0:
            readable_number = f"{number}bp"
        else:
            giga_num = int(number / 1e9)
            readable_number = f"{giga_num}Gbp"
    else:
        readable_number = f"{number}bp"

    return readable_number


def prepare_summary(args, stats, proc_timings):

    proc_timings = np.sort(np.array(proc_timings, dtype=float))
    proc_median = round(proc_timings[proc_timings.size//2], 6)
    proc_mean = round(proc_timings.mean(), 6)
    proc_top = round(proc_timings[int(proc_timings.size * 0.99)], 6)

    ref_size = args.coverage_ref_size
    length_thresholds = args.length_thresholds
    if min(length_thresholds) != 0:
        length_thresholds = [0] + length_thresholds

    summary_stats = [
        ("all", "sec_per_seq_median", proc_median),
        ("all", "sec_per_seq_mean", proc_mean),
        ("all", "sec_per_seq_99pct", proc_top)
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


def write_output_file(file_path, data, dump_index):

    if file_path is not None:
        file_path.parent.mkdir(exist_ok=True, parents=True)
        data.to_csv(
            file_path, sep="\t",
            header=True, index=dump_index
        )
    return


def main():

    args = parse_command_line()

    all_inputs = args.input_files
    if isinstance(all_inputs, list) and len(all_inputs) == 1:
        if all_inputs[0].suffix in ["fofn", "fof", "fofp"]:
            all_inputs = []
            with open(all_inputs[0], "r") as fofn:
                for line in fofn:
                    this_file = pl.Path(line.strip()).resolve(strict=True)
                    all_inputs.append(this_file)

    processors = setup_sequence_processors(args)

    proc_sequence = functools.partial(
        process_sequence,
        processors
    )

    stats = []
    index_records = []
    proc_timings = []
    with mp.Pool(args.cores) as pool:
        for input_file in sorted(all_inputs):
            file_name = input_file.name
            with dnaio.open(str(input_file)) as seq_file:
                stats_iter = pool.imap_unordered(proc_sequence, seq_file)
                for seq_name, seq_hash, seq_stats, proc_time in stats_iter:
                    index_records.append((seq_name, file_name, seq_hash))
                    stats.append(seq_stats)
                    proc_timings.append(proc_time)

    proc_timings = pd.DataFrame.from_records(proc_timings, index="seq_name")
    write_output_file(args.output_timings, proc_timings, True)

    mindex = pd.MultiIndex.from_tuples(index_records, names=["seq_name", "seq_source", "seq_hash"])
    stats = pd.DataFrame.from_records(stats, index=mindex)
    stats.fillna(0, inplace=True)
    write_output_file(args.output_statistics, stats, True)

    summary = prepare_summary(args, stats, proc_timings["total"].values)
    write_output_file(args.output_summary, summary, False)

    return 0



if __name__ == "__main__":
    main()
