import pathlib
import hashlib
import collections
import re

import pandas

SAMPLES = None
MAP_SAMPLE_TO_INPUT_FILES = None

TRIO_SAMPLES = None
UNPHASED_SAMPLES = None

CONSTRAINT_TRIO_SAMPLES = None
CONSTRAINT_UNPHASED_SAMPLES = None


def process_sample_sheet():

    SAMPLE_SHEET_FILE = pathlib.Path(config["samples"]).resolve(strict=True)
    SAMPLE_SHEET = pandas.read_csv(
        SAMPLE_SHEET_FILE,
        sep="\t",
        header=0
    )

    # step 1: each row is a sample,
    # just collect the input files
    sample_input, trio_samples, unphased_samples = collect_input_files(SAMPLE_SHEET)
    all_samples = sorted(sample_input.keys())
    assert len(all_samples) == SAMPLE_SHEET.shape[0]

    global SAMPLES
    SAMPLES = all_samples
    global TRIO_SAMPLES
    TRIO_SAMPLES = trio_samples
    global UNPHASED_SAMPLES
    UNPHASED_SAMPLES = unphased_samples

    global MAP_SAMPLE_TO_INPUT_FILES
    MAP_SAMPLE_TO_INPUT_FILES = sample_input

    global CONSTRAINT_TRIO_SAMPLES
    CONSTRAINT_TRIO_SAMPLES = _build_constraint(trio_samples)
    global CONSTRAINT_UNPHASED_SAMPLES
    CONSTRAINT_UNPHASED_SAMPLES = _build_constraint(unphased_samples)

    return


def collect_input_files(sample_sheet):
    """
    The output of this function should
    be sufficient to run the workflow
    in single-sample mode
    """
    sample_input = collections.defaultdict(dict)
    trio_samples = set()
    unphased_samples = set()

    for row in sample_sheet.itertuples():
        sample = row.sample
        hifi_input, hifi_hashes = collect_sequence_input(row.hifi)
        ont_input, ont_hashes = collect_sequence_input(row.ont)
        if row.target == "trio":
            trio_samples.add(sample)
            hap1_db = row.hap1
            assert hap1_db.endswith("meryl")
            assert pathlib.Path(hap1_db).resolve(strict=True).is_dir()
            sample_input[sample]["hap1"] = hap1_db

            hap2_db = row.hap2
            assert hap2_db.endswith("meryl")
            assert pathlib.Path(hap2_db).resolve(strict=True).is_dir()
            sample_input[sample]["hap2"] = hap2_db

            assert sample_input[sample]["hap1"] != sample_input[sample]["hap2"]

        elif row.target == "unphased":
            unphased_samples.add(sample)
        else:
            raise ValueError(row.target)
        sample_input[sample]["hifi"] = hifi_input
        sample_input[sample]["ont"] = ont_input

    return sample_input, trio_samples, unphased_samples


def _read_input_files_from_fofn(fofn_path):
    """Read input file listing from
    file of file names
    TODO: candidate for inclusion in template
    """

    input_files = []
    with open(fofn_path, "w") as listing:
        for line in listing:
            if not line.strip():
                continue
            try:
                file_path = pathlib.Path(line.strip()).resolve(strict=True)
            except FileNotFoundError:
                try:
                    file_path = DATA_ROOT.joinpath(line.strip()).resolve(strict=True)
                except FileNotFoundError:
                    err_msg = "\nERROR\n"
                    err_msg += f"Cannot find file: {line.strip}\n"
                    err_msg += f"Data root is set to: {DATA_ROOT}\n"
                    sys.stderr.write(err_msg)
                    raise
            input_files.append(file_path)

    return sorted(input_files)


def subset_path(full_path):
    """This helper exists to reduce
    the absolute path to a file
    to just the file name and its
    parent.
    TODO: should be codified as part
    of the template utilities to improve
    infrastructure portability of active
    workflows
    """
    folder_name = full_path.parent.name
    file_name = full_path.name
    subset_path = f"{folder_name}/{file_name}"
    # if it so happens that the file resides
    # in a root-level location, strip off
    # leading slash
    return subset_path.strip("/")


def collect_sequence_input(path_spec):
    """
    Generic function to collect HiFi or ONT/Nanopore
    input (read) files
    """
    input_files = []
    input_hashes = []
    for sub_input in path_spec.split(","):
        input_path = pathlib.Path(sub_input).resolve()
        if input_path.is_file() and input_path.name.endswith(".fofn"):
            fofn_files = _read_input_files_from_fofn(input_path)
            fofn_files = [
                hashlib.sha256(
                    subset_path(fp).encode("utf-8")
                ).hexdigest() for fp in fofn_files
            ]
            input_files.extend(collected_files)
            input_hashes.extend(collected_hashes)
        elif input_path.is_file():
            input_hash = hashlib.sha256(
                subset_path(input_path).encode("utf-8")
            ).hexdigest()
            input_files.append(input_path)
            input_hashes.append(input_hash)
        elif input_path.is_dir():
            collected_files = _collect_files(input_path)
            collected_hashes = [
                hashlib.sha256(
                    subset_path(fp).encode("utf-8")
                ).hexdigest() for fp in collected_files
            ]
            input_files.extend(collected_files)
            input_hashes.extend(collected_hashes)
        else:
            raise ValueError(f"Cannot handle input: {sub_input}")
    return input_files, input_hashes


def _build_constraint(values):
    escaped_values = sorted(map(re.escape, map(str, values)))
    constraint = "(" + "|".join(escaped_values) + ")"
    return constraint


def _collect_files(folder):

    all_files = set()
    for pattern in config["input_file_ext"]:
        pattern_files = set(folder.glob(f"**/*.{pattern}"))
        all_files = all_files.union(pattern_files)
    all_files = [f for f in sorted(all_files) if f.is_file()]
    if len(all_files) < 1:
        raise ValueError(f"No input files found underneath {folder}")
    return all_files

process_sample_sheet()
