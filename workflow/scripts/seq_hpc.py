#!/usr/bin/env python3

import argparse as argp
import contextlib as ctl
import io
import pathlib as pl
import sys
import time

import dnaio
import xopen


def parse_command_line():

    parser = argp.ArgumentParser(prog="seq_hpc.py")

    parser.add_argument(
        "--input",
        "-in",
        "-i",
        default="stdin",
        type=str,
        dest="input",
        help="Full path to FASTA/FASTQ input file or 'stdin'. Default: stdin"
    )

    parser.add_argument(
        "--in-format",
        "-fmt",
        default="fasta",
        type=str,
        choices=["fasta", "fastq"],
        dest="input_format",
    )

    parser.add_argument(
        "--force-format",
        "-f",
        action="store_true",
        default=False,
        dest="force_format",
        help="Only for file input: force format, otherwise use auto-detection."
    )

    parser.add_argument(
        "--output",
        "-out",
        "-o",
        default="stdout",
        type=str,
        dest="output",
        help="Full path to FASTA output file or 'stdout'. Default: stdout"
    )

    parser.add_argument(
        "--report",
        "-r",
        action="store_true",
        default=False,
        dest="report",
        help="If set, write brief report to sys.stderr. Default: False"
    )

    parser.add_argument(
        "--buffer-size",
        "-b",
        type=int,
        default=0,
        dest="buffer_size",
        help="Set buffer size (#char) or 0 for no buffer. Default: 0"
    )

    args = parser.parse_args()

    return args


def homopolymer_compress(sequence):
    """NB: this would fail on empty input

    Args:
        sequence (str): sequence to hpc
    Returns:
        str: hpc sequence
    """
    uncompressed_length = len(sequence)
    hpc_seq = "".join(a if a != b else "" for a,b in zip(sequence[:-1], sequence[1:]))
    if not hpc_seq:
        hpc_seq = sequence[0]
    if hpc_seq[-1] != sequence[-1]:
        hpc_seq += sequence[-1]
    compressed_length = len(hpc_seq)
    ratio = compressed_length / uncompressed_length
    assert ratio <= 1
    return hpc_seq, ratio


def check_implementation():

    test_strings = [
        "AAA", "ACGT", "ACAA", "XXXNNX",
        "AACGT", "ACCCCGGGGT","ACGTTT"
    ]

    result_strings = [
        "A", "ACGT", "ACA", "XNX",
        "ACGT", "ACGT", "ACGT"
    ]

    assert all(
        [homopolymer_compress(test) == result]
        for test, result in zip(test_strings, result_strings)
    )

    return


def main():

    _ = check_implementation()

    args = parse_command_line()

    if args.input not in ["stdin", "-", "/dev/stdin", ""]:
        infile = pl.Path(args.input).resolve(strict=True)
    else:
        infile = sys.stdin.buffer

    if not args.force_format:
        input_format = None
    else:
        input_format = args.input_format

    if args.output not in ["stdout", "-", "/dev/stdout", ""]:
        outfile = pl.Path(args.output).resolve()
        outfile.parent.mkdir(exist_ok=True, parents=True)
    else:
        outfile = sys.stdout.buffer


    buffer_limit = args.buffer_size
    buffered = 0
    use_buffer = buffer_limit > 0
    buffer_obj = None
    write_buffer = None

    avg_hpc_ratio = 0.
    process_start = time.perf_counter()
    with ctl.ExitStack() as exs:
        read_input = exs.enter_context(
            dnaio.open(infile, fileformat=input_format, mode="r")
        )

        if use_buffer:
            buffer_obj = io.BytesIO()
            write_buffer = dnaio.open(buffer_obj, mode="w", fileformat="fasta")
            write_output = exs.enter_context(
                xopen.xopen(outfile, mode="wb", compresslevel=5)
            )
        else:
            write_output = exs.enter_context(
                dnaio.open(
                    outfile, fileformat="fasta", mode="w",
                    qualities=False, compression_level=5
                )
            )

        hpc = homopolymer_compress
        for record_num, seq_record in enumerate(read_input, start=1):
            hpc_seq, ratio = hpc(seq_record.sequence.upper())
            avg_hpc_ratio += ratio

            if use_buffer:
                write_buffer.write(seq_record.name, hpc_seq)
                buffered += len(seq_record.name)
                buffered += len(hpc_seq)
                if buffered > buffer_limit:
                    write_output.write(buffer_obj.getvalue())
                    buffer_obj = io.BytesIO()
                    write_buffer = dnaio.open(buffer_obj, mode="w", fileformat="fasta")
                    buffered = 0
            else:
                write_output.write(seq_record.name, hpc_seq)

        if buffered > 0:
            write_output.write(buffer_obj.getvalue())

    process_end = time.perf_counter()
    total_time = round(process_end - process_start, 3)

    avg_hpc_ratio = round(avg_hpc_ratio / record_num, 3)
    if args.report:
        sys.stderr.write(
            f"\n=== seq_hpc.py report ==="
            f"\nProcessed records: {record_num}"
            f"\nAverage homopolymer-compression ratio: {avg_hpc_ratio}"
            f"\nTotal processing time: ~{total_time} sec\n\n"
        )

    return 0


if __name__ == "__main__":
    main()
