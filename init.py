#!/usr/bin/env python3

import argparse as argp
import logging
import logging.config as logconf
import os
import pathlib as pl
import subprocess as sp
import sys


def create_execution_environment(repo_folder, project_folder, dev_only):
    """Create Conda environments if any Conda-like executable
    is found on $PATH

    Args:
        repo_folder (pathlib.Path): This repository checkout location
        project_folder (pathlib.Path): Project folder, assumed to be
        one above the repo folder
        dev_only (boolean): Only create Conda environment for
        development purposes

    Raises:
        RuntimeError: In dev only mode, a Conda-like executable
        must be available

        subprocess.CalledProcessError: Propagated if Conda
        environment cannot be created

    Returns:
        None: placeholder
    """

    logger = logging.getLogger(__name__)
    check_executables = ["mamba", "conda"]
    use_executable = None
    for executable in check_executables:
        try:
            _ = sp.check_call(
                [executable, "--version"],
                shell=False,
                stdout=sp.DEVNULL,
                stderr=sp.DEVNULL,
            )
            use_executable = executable
            break
        except sp.CalledProcessError as spe:
            logger.warning(f"Executable {executable} not available: {spe}")
    if use_executable is None:
        logger.warning(
            "No executable available to create execution (conda) environment"
        )
        if dev_only:
            raise RuntimeError("No Conda executable available, cannot create dev env")
    else:
        logger.debug(
            f"Found Conda executable {use_executable} - creating environment..."
        )
        if dev_only:
            logger.debug('Development mode set, select "dev_env.yaml" file.')
            yaml_file = repo_folder / pl.Path("workflow", "envs", "dev_env.yaml")
            yaml_file = yaml_file.resolve(strict=True)
            env_prefix = repo_folder / pl.Path("dev_env")
        else:
            yaml_file = repo_folder / pl.Path("workflow", "envs", "exec_env.yaml")
            yaml_file = yaml_file.resolve(strict=True)
            env_prefix = project_folder / pl.Path("exec_env")
        logger.info(
            f"Creating Snakemake execution environment at location: {env_prefix}"
        )
        logger.debug("Setting up the execution environment may take a while...")
        call_args = [
            use_executable,
            "env",
            "create",
            "--quiet",
            "--force",
            "-f",
            str(yaml_file),
            "-p",
            str(env_prefix),
        ]
        try:
            proc_out = sp.run(call_args, shell=False, capture_output=True, check=False)
            proc_out.check_returncode()  # check after to get stdout/stderr
        except sp.CalledProcessError as spe:
            logger.error(f"Could not create Snakemake execution environment: {spe}")
            logger.error(f"\n=== STDOUT ===\n{proc_out.stdout.decode('utf-8')}")
            logger.error(f"\n=== STDERR ===\n{proc_out.stderr.decode('utf-8')}")
            raise
    return None


def setup_logging(project_dir, debug_mode, dev_only):
    """Setup logging to stderr stream and file

    Args:
        project_dir (pathlib.Path): Project folder,
        assumed to be one above repo location
        debug_mode (boolean): log verbose / debug
        dev_only (boolean): create only Conda dev
        environment, do not create init.log file

    Returns:
        pathlib.Path or os.devnull: log file location
    """

    base_level = "DEBUG" if debug_mode else "WARNING"
    log_file_location = project_dir / pl.Path("init.log")
    if dev_only:
        log_file_location = os.devnull

    log_config = {
        "version": 1,
        "root": {"handlers": ["stream", "file"], "level": "DEBUG"},
        "handlers": {
            "stream": {
                "formatter": "default",
                "class": "logging.StreamHandler",
                "level": base_level,
                "stream": sys.stderr,
            },
            "file": {
                "formatter": "default",
                "class": "logging.FileHandler",
                "level": "INFO",
                "filename": log_file_location,
            },
        },
        "formatters": {
            "default": {
                "format": "%(asctime)s : \
                    %(levelname)s - \
                    %(funcName)s - \
                    ln:%(lineno)d >> \
                    %(message)s",
                "datefmt": "%Y-%m-%d %H:%M:%S",
            },
        },
    }
    logconf.dictConfig(log_config)
    return log_file_location


def parse_command_line():
    """Create command line parser

    Returns:
        argparse.Namespace: command line options
    """
    parser = argp.ArgumentParser()
    parser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Print log messages to stderr.",
        dest="debug",
    )
    parser.add_argument(
        "--dev-only",
        action="store_true",
        default=False,
        help="Only create a Conda environment for development purposes \
            (no working directory hierarchy).",
        dest="dev_only",
    )
    args = parser.parse_args()
    return args


def create_wd_folders(project_dir):
    """Create folder hierarchy starting
    at Snakemake's future working directory

    Args:
        project_dir (pathlib.Path): Project folder,
        assumed to be one above repo location

    Returns:
        None: placeholder
    """

    logger = logging.getLogger(__name__)
    logger.info("Creating Snakemake working directory structure")
    wd_toplevel = project_dir / pl.Path("wd")
    wd_toplevel.mkdir(exist_ok=True, parents=True)

    subfolders = [
        ("proc",),
        ("results",),
        ("log",),
        ("rsrc",),
        ("log", "cluster_jobs", "err"),
        ("log", "cluster_jobs", "out"),
        ("global_ref",),
        ("local_ref",),
    ]

    for sub in subfolders:
        full_path = wd_toplevel / pl.Path(*sub)
        logger.info(f"Creating path {full_path}")
        full_path.mkdir(exist_ok=True, parents=True)

    return None


def main():
    """Main function

    Returns:
        integer: explicit 0 on success
    """
    args = parse_command_line()
    repo_location = pl.Path(__file__).resolve(strict=True).parent
    project_dir = repo_location.parent
    log_file_location = setup_logging(project_dir, args.debug, args.dev_only)
    logger = logging.getLogger(__name__)
    logger.info(f"Repository location: {repo_location}")
    logger.info(f"Project directory: {project_dir}")
    logger.info(f"Log file location: {log_file_location}")
    create_execution_environment(repo_location, project_dir, args.dev_only)
    if not args.dev_only:
        create_wd_folders(project_dir)

    return 0


if __name__ == "__main__":
    main()
