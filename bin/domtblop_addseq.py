#!/usr/bin/env python

"""
Add sequences to each Query Result that has maching ids in a fasta file.

Usage:
    domtblop.py group [options] <serialized_domtblout> <fasta_file>


Arguments:
    <serialized_domtblout> STR      Path to the serialized domtblout file.
                                    Use "-" to read from stdin. (default: "-")
    <fasta_file> STR                Path to the fasta file containing the sequences.


Options:
    --nucleotide                    Include sequence as nucleotide sequence.
    -h, --help                      Show this help message and exit
    -l, --loglevel STR              Set the logging level [default: INFO]
"""

import argparse
import json
import logging
import os
import sys
from typing import (
        Dict,
        List,
        Tuple,
        )

from Bio.SeqIO.FastaIO import SimpleFastaParser

from domtblop_parser import (
        HmmscanQueryResult,
        UnexpectedQueryIdFormat,
        CustomEncoder,
        )
from domtblop_utils import (
        setup_logger,
        read_input,
    )


def setup_argparse() -> argparse.ArgumentParser:
    """
    Sets up the argparse instance for command-line arguments.

    Returns:
    argparse.ArgumentParser: Configured ArgumentParser instance.
    """

    parser = argparse.ArgumentParser(
        add_help=False,
    )

    # Required arguments
    parser.add_argument(
        "serialized_domtblout",
        metavar="<serialized_domtblout>",
        type=str,
        nargs="?",
        default="-",
        )
    parser.add_argument(
        "fasta_file",
        metavar="<fasta_file>",
        type=str,
        nargs=1,
        )

    # Optional arguments
    parser.add_argument(
        "--nucleotide",
        action="store_true",
        default=False,
        help="Include sequence as parent sequence."
    )
    parser.add_argument(
        "-l", "--loglevel",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )
    parser.add_argument(
        "-h", "--help",
        action="store_true",
        default=False
    )

    return parser


def setup_config(args: List[str],) -> Tuple[
        argparse.Namespace,
        logging.Logger
        ]:
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

    if not os.path.exists(config.serialized_domtblout) and config.serialized_domtblout != "-":
        logger.error(
            f"Serialized domtblout file '{config.serialized_domtblout}' "
            "does not exist."
        )
        raise FileNotFoundError

    if not os.path.exists(config.fasta_file[0]):
        logger.error(
            f"FASTA file '{config.fasta_file[0]}' does not exist."
        )
        raise FileNotFoundError


    return config, logger


def parse_fasta(fasta_file: str) -> Dict[str, str]:
    """
    Parse a fasta file and return a dictionary with the title as key and the
    sequence as value.

    Args:
        fasta_file (str): Path to the fasta file.

    Returns:
        Dict[str, str]: Dictionary with the title as key and the sequence as value.

    Raises:
        FileNotFoundError: If the fasta file does not exist.
    """

    fasta_contents = {}
    with open(fasta_file) as handle:
        for title, seq in SimpleFastaParser(handle):
            fasta_contents[title] = seq

    return fasta_contents


def run(args: List[str]) -> None:

    try:
        config, logger = setup_config(args)
    except (ValueError, FileNotFoundError):
        sys.exit(1)
    logger.info(
        "Running with the following configuration: " 
        + ", ".join(f"{k}={v}" for k,v in config.__dict__.items())
        )


    try:
        fasta_contents = parse_fasta(config.fasta_file[0])
    except FileNotFoundError as e:
        logger.error(f"{e}")
        sys.exit(1)

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

        if config.nucleotide:
            query_result.nucleotide_sequence = fasta_contents.get(query_result.query_id, None)

        else:
            query_result.aminoacid_sequence = fasta_contents.get(query_result.query_id, None)

            query_result.mk_domain_alignment_fragment_sequences()

        try:
            print(json.dumps(query_result, cls=CustomEncoder))

        except BrokenPipeError:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())
            sys.exit(1)

