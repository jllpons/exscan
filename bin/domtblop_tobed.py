#!/usr/bin/env python

"""
Write structured query results to stdout in BED format.

Usage:
    domtblop.py tobed [options] <serialized_domtblout>


Arguments:
    <serialized_domtblout> STR      Path to the serialized domtblout file.
                                    Use "-" to read from stdin. (default: "-")

Options:
    -h, --help                      Show this help message and exit
    -l, --loglevel STR              Set the logging level [default: INFO]
"""

import argparse
from dataclasses import dataclass
import json
import logging
import os
import sys
from typing import List, Tuple

from domtblop_parser import (
    HmmscanQueryResult,
    UnexpectedQueryIdFormat,
)
from domtblop_utils import setup_logger, read_input


@dataclass
class BedFeature:
    seqid: str
    start: int
    end: int
    name: str

    def __str__(self):
        return f"{self.seqid}\t{self.start}\t{self.end}\t{self.name}"


def bed_from_query_result(query_result: HmmscanQueryResult) -> BedFeature:
    """
    Convert a HmmscanQueryResult object to a list of BED entries.

    Args:
        query_result (HmmscanQueryResult): A HmmscanQueryResult object.

    Returns:
        BedFeature: A BedFeature object.
    """
    seqid = query_result.parent_sequence.sequence_id
    start = (
        query_result.parent_sequence.start - 1
        if query_result.parent_sequence.start > 0
        else 0
    )
    end = query_result.parent_sequence.end
    query_id = query_result.query_id

    feature = BedFeature(seqid=seqid, start=start, end=end, name=query_id)

    return feature


def setup_argparse() -> argparse.ArgumentParser:
    """
    Sets up the argparse instance for command-line arguments.

    Returns:
        argparse.ArgumentParser: Configured ArgumentParser instance.
    """

    parser = argparse.ArgumentParser(add_help=False)

    # Required arguments
    parser.add_argument(
        "serialized_domtblout",
        metavar="<serialized_domtblout>",
        type=str,
        nargs="?",
        default="-",
    )

    # Optional arguments
    parser.add_argument(
        "--out-format",
        type=str,
        default="bed",
        choices=["bed", "gff3"],
    )
    parser.add_argument(
        "-l",
        "--loglevel",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )
    parser.add_argument("-h", "--help", action="store_true", default=False)

    return parser


def setup_config(
    args: List[str],
) -> Tuple[argparse.Namespace, logging.Logger]:
    """
    Setup configuration for the script.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = setup_argparse()
    config = parser.parse_args(args)

    if len(sys.argv) == 1 and sys.stdin.isatty():
        print(__doc__)
        sys.exit(1)

    if config.help:
        print(__doc__)
        sys.exit(0)

    logger = setup_logger(config.loglevel)

    if not config.serialized_domtblout:
        logger.error("No serialized domtblout file provided.")
        raise ValueError

    if (
        not os.path.exists(config.serialized_domtblout)
        and config.serialized_domtblout != "-"
    ):
        logger.error(
            f"Serialized domtblout file '{config.serialized_domtblout}' "
            "does not exist."
        )
        raise FileNotFoundError

    return config, logger


def run(args: List[str]) -> None:
    try:
        config, logger = setup_config(args)
    except (ValueError, FileNotFoundError):
        sys.exit(1)
    logger.info(
        "Running with the following configuration: "
        + ", ".join(f"{k}={v}" for k, v in config.__dict__.items())
    )

    try:
        file_handle = read_input(config.serialized_domtblout)
    except FileNotFoundError:
        sys.exit(1)

    for line in file_handle:
        try:
            query_result = HmmscanQueryResult.from_json(json.loads(line))
        except UnexpectedQueryIdFormat as e:
            logger.error(e)
            sys.exit(1)

        feature = bed_from_query_result(query_result)

        try:
            print(feature)
        except BrokenPipeError:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())
            sys.exit(1)


if __name__ == "__main__":
    run(sys.argv[1:])
