import pathlib
import hashlib
import collections
import re

import pandas

SAMPLES = None
MAP_SAMPLE_TO_INPUT_FILE = None
MAP_SAMPLE_TO_PART_IDS = None

CONSTRAINT_SAMPLES_SINGLE_INPUT = None
CONSTRAINT_SAMPLES_MULTI_INPUT = None

CHILDREN = None
MAP_TRIOS = None


def process_sample_sheet():

    SAMPLE_SHEET_FILE = pathlib.Path(config["samples"]).resolve(strict=True)
    SAMPLE_SHEET = pandas.read_csv(
        SAMPLE_SHEET_FILE,
        sep="\t",
        header=0
    )

    # step 1: each row is a sample,
    # just collect the input files
    sample_input, single_input_samples, multi_input_samples = collect_input_files(SAMPLE_SHEET)
    all_samples = single_input_samples.union(multi_input_samples)
    assert len(all_samples) == SAMPLE_SHEET.shape[0]

    global SAMPLES
    SAMPLES = sorted(all_samples)
    global MAP_SAMPLE_TO_INPUT_FILE
    MAP_SAMPLE_TO_INPUT_FILE = sample_input
    global MAP_SAMPLE_TO_PART_IDS
    MAP_SAMPLE_TO_PART_IDS = dict(
        (sample, path_hash) for (sample, path_hash) in sample_input.keys() if path_hash is not None
    )
    global CONSTRAINT_SAMPLES_SINGLE_INPUT
    CONSTRAINT_SAMPLES_SINGLE_INPUT = _build_constraint(single_input_samples)
    global CONSTRAINT_SAMPLES_MULTI_INPUT
    CONSTRAINT_SAMPLES_MULTI_INPUT = _build_constraint(multi_input_samples)

    # step 2: check if any families/trios
    # are part of the sample sheet
    children, trios = collect_family_infos(SAMPLE_SHEET)

    global CHILDREN
    CHILDREN = children
    global MAP_TRIOS
    MAP_TRIOS = trios

    return


def collect_input_files(sample_sheet):
    """
    The output of this function should
    be sufficient to run the workflow
    in single-sample mode (i.e. ignoring
    trio information)
    """
    sample_input = dict()
    single_input_samples = set()
    multi_input_samples = set()

    for row in sample_sheet.itertuples():
        input_files = []
        input_hashes = []
        for sub_input in row.input.split(","):
            input_path = pathlib.Path(sub_input).resolve()
            if input_path.is_file():
                input_hash = hashlib.sha256(str(input_path).encode("utf-8")).hexdigest()
                input_files.append(input_path)
                input_hashes.append(input_hash)
            elif input_path.is_dir():
                collected_files = _collect_files(input_path)
                collected_hashes = [
                    hashlib.sha256(str(f).encode("utf-8")).hexdigest() for f in collected_files
                ]
                input_files.extend(collected_files)
                input_hashes.extend(collected_hashes)
            else:
                raise ValueError(f"Cannot handle input: {sub_input}")


        if len(input_files) > 1:
            for infile, path_hash in zip(input_files, input_hashes):
                sample_input[(row.sample, path_hash)] = infile
            multi_input_samples.add(row.sample)
        else:
            sample_input[(row.sample, None)] = input_files[0]
            single_input_samples.add(row.sample)
    return sample_input, single_input_samples, multi_input_samples


def collect_family_infos(sample_sheet):

    children = set()
    families = collections.defaultdict(dict)
    proper_trios = collections.defaultdict(dict)

    for row in sample_sheet.itertuples():
        if not hasattr(row, "family") or pandas.isnull(row.family):
            continue
        families[row.family][row.member] = row.sample

    for family, members in families.items():
        if len(members) != 3:
            raise ValueError(f"Family must have 3 members: {members}")
        child = members["child"]
        if child in proper_trios:
            raise ValueError(f"Child member of several trios: {child}")
        proper_trios[child]["mother"] = members["mother"]
        proper_trios[child]["father"] = members["father"]
        children.add(child)

    return children, proper_trios


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
