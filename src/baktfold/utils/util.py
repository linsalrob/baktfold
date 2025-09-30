import os
import shutil
import sys
import time
from pathlib import Path
from typing import Any, Dict

import click
from loguru import logger


class OrderedCommands(click.Group):
    """This class will preserve the order of subcommands, which is useful when printing --help"""

    def list_commands(self, ctx: click.Context):
        return list(self.commands)


def baktfold_base(rel_path):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    with open(baktfold_base("VERSION"), "r") as f:
        version = f.readline()
    return version


def echo_click(msg, log=None):
    click.echo(msg, nl=False, err=True)
    if log:
        with open(log, "a") as lo:
            lo.write(msg)


def print_citation():
    with open(baktfold_base("CITATION"), "r") as f:
        for line in f:
            echo_click(line)


log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)

"""
begin and end functions
"""


def begin_baktfold(params: Dict[str, Any], subcommand: str) -> int:
    """
    Begin baktfold process.

    Parameters:
        params (Dict[str, Any]): A dictionary of parameters for baktfold.
        subcommand (str): Subcommand indicating the baktfold operation.

    Returns:
        int: Start time of the baktfold process.
    """
    # get start time
    start_time = time.time()
    # initial logging stuff
    log_file = os.path.join(params["--output"], f"baktfold_{subcommand}_{start_time}.log")
    # adds log file
    logger.add(log_file)
    logger.add(lambda _: sys.exit(1), level="ERROR")

    print_splash()
    logger.info("baktfold: rapid & standardized annotation of bacterial genomes, MAGs & plasmids using protein structural information")

    logger.info(f"You are using baktfold version {get_version()}")
    logger.info("Repository homepage is https://github.com/gbouras13/baktfold")
    logger.info(f"You are running baktfold {subcommand}")
    logger.info(f"Listing parameters")
    for key, value in params.items():
        logger.info(f"Parameter: {key} {value}")

    return start_time


def end_baktfold(start_time: float, subcommand: str) -> None:
    """
    Finish baktfold process and log elapsed time.

    Parameters:
        start_time (float): Start time of the process.
        subcommand (str): Subcommand name indicating the baktfold operation.

    Returns:
        None
    """

    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info(f"baktfold {subcommand} has finished")
    logger.info("Elapsed time: " + str(elapsed_time) + " seconds")


# need the logo here eventually
def print_splash():
    click.echo(
        """\b

  _           _    _    __      _     _ 
 | |         | |  | |  / _|    | |   | |
 | |__   __ _| | _| |_| |_ ___ | | __| |
 | '_ \ / _` | |/ / __|  _/ _ \| |/ _` |
 | |_) | (_| |   <| |_| || (_) | | (_| |
 |_.__/ \__,_|_|\_\\__|_| \___/|_|\__,_|
                                        
                                        
"""
    )


def remove_file(file_path: Path) -> None:
    """
    Remove a file if it exists.

    Parameters:
        file_path (Path): Path to the file to remove.

    Returns:
        None
    """
    if file_path.exists():
        file_path.unlink()  # Use unlink to remove the file


def remove_directory(dir_path: Path) -> None:
    """
    Remove a directory and all its contents if it exists.

    Parameters:
        dir_path (Path): Path to the directory to remove.

    Returns:
        None
    """
    if dir_path.exists():
        shutil.rmtree(dir_path)


def touch_file(path: Path) -> None:
    """
    Update the access and modification times of a file to the current time, creating the file if it does not exist.

    Parameters:
        path (Path): Path to the file.

    Returns:
        None
    """
    with open(path, "a"):
        os.utime(path, None)


def clean_up_temporary_files(output: Path) -> None:
    """
    Clean up temporary files generated during the baktfold process.

    Parameters:
        output (Path): Path to the output directory.

    Returns:
        None
    """
    result_high_tsv: Path = Path(output) / "foldseek_results_high.tsv"
    result_low_tsv: Path = Path(output) / "foldseek_results_low.tsv"
    result_tsv: Path = Path(output) / "foldseek_results.tsv"
    result_tsv_custom: Path = Path(output) / "foldseek_results_custom.tsv"
    foldseek_db: Path = Path(output) / "foldseek_db"
    result_db_base: Path = Path(output) / "result_db"
    temp_db: Path = Path(output) / "temp_db"
    aln_db: Path = Path(output) / "aln_db"
    remove_directory(result_db_base)
    remove_directory(temp_db)
    remove_directory(foldseek_db)
    remove_directory(aln_db)
    remove_file(result_tsv)
    remove_file(result_tsv_custom)
    remove_file(result_high_tsv)
    remove_file(result_low_tsv)
