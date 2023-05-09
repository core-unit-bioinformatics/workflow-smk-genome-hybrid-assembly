
def assemble_verkko_screen_string():
    """This helper function exists to
    properly assemble the command line
    string for Verkko's >screen< option,
    which only has a default for human.
    """

    screen_files_exist = config.get("verkko_screen_files_exist", True)

    screen_opt = ""
    screen_spec = config.get("verkko_screen", "")
    if screen_spec == "human":
        screen_opt = " --screen human "
    elif isinstance(screen_spec, list):
        assert isinstance(screen_spec[0], dict)
        screen_opt = " "
        for spec in screen_spec:
            label, fasta = [(k, v) for k, v in spec.items()][0]
            fasta = pathlib.Path(fasta).resolve(strict=screen_files_exist)
            screen_opt += f"--screen {label} {fasta} "
    elif isinstance(screen_spec, dict):
        screen_opt = " "
        for label, fasta in screen_spec.items():
            fasta = pathlib.Path(fasta).resolve(strict=screen_files_exist)
            screen_opt += f"--screen {label} {fasta} "
    else:
        raise ValueError(f"Cannot parse Verkko screen specification: {screen_spec}")

    if VERBOSE:
        logout(f"Verkko screen option set to: {screen_opt}")

    return screen_opt


def increase_mbg_resources(attempt):
    """For some HiFi datasets, MBG requires
    much more time to build the initial graph.
    This helper function exists to directly
    increase the MBG resource requirements
    if the Verkko run is restarted.

    In Verkko v1.3.1, this pertains to ...
    ... option string `--mbg-run`
    and affects settings ...
    #  build-graph
    mbg_n_cpus:          '4'
    mbg_mem_gb:          '0'
    mbg_time_h:          '72'
    ... for rule `1-buildGraph.sm::buildGraph`
    """
    mbg_resources = ""
    if int(attempt) > 1:
        # this is CPU - MEM_GB - TIME_HRS
        mbg_resources = "--mbg-run 8 160 72"
    return mbg_resources


def increase_process_ont_resources(attempt):
    """Singular cases where processing the ONT
    paths is OOM killed. Use the same strategy
    as above and increase the initial job
    resources if Verkko is restarted
    (NB: Verkko also restarts its jobs)

    In Verkko v1.3.1, this pertains to...
    ... option string `--pop-run`
    and affects settings ...
    #  process_ont_paths
    pop_n_cpus:          '1'
    pop_mem_gb:          '64'
    pop_time_h:          '24'
    ... for rule `4-processONT.sm::processONT`
    """
    pop_resources = ""
    if int(attempt) > 1:
        # this is CPU - MEM_GB - TIME_HRS
        pop_resources = "--pop-run 1 96 48"
    return pop_resources

