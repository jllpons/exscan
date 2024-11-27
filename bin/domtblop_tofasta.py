#!/usr/bin/env python

"""
Write structured query results to stdout in FASTA format.

Usage:
    domtblop.py tofasta [options] <serialized_domtblout>


Arguments:
    <serialized_domtblout> STR      Path to the serialized domtblout file.
                                    Use "-" to read from stdin. (default: "-")

    --seq-kind STR                  Kind of sequence to extract from the serialized domtblout.
                                    Choices:
                                        - "nucleotide": The nucleotide sequence of the ORFs that contains the domain hit.
                                        - "protein": The protein sequence of the ORFs that contains the domain hit.
                                        - "domain-alignment": The protein sequence of the ORFs that could be aligned against the domain.

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
class Fasta:
    header: str
    sequence: str

    def __str__(self):
        seq_newline = "\n".join(self.sequence[i:i+80] for i in range(0, len(self.sequence), 80))
        return f">{self.header}\n{seq_newline}"


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
    parser.add_argument(
        "--seq-kind",
        type=str,
        choices=["nucleotide", "protein", "domain-alignment"],
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


def setup_config(args: List[str],) -> Tuple[argparse.Namespace, logging.Logger]:
    """
    Setup configuration for the script.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = setup_argparse()
    config = parser.parse_args(args)

    if len(sys.argv) == 2 and sys.stdin.isatty():
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
            f"Serialized domtblout file '{config.serialized_domtblout}' " "does not exist."
        )
        raise FileNotFoundError


    if not config.seq_kind:
        logger.error("No sequence kind provided. Run with -h/--help for more information.")
        raise ValueError


    if config.seq_kind not in ["nucleotide", "protein", "domain-alignment"]:
        logger.error(
            f"Sequence kind '{config.seq_kind}' not recognized. "
            "Choose from {'nucleotide', 'protein', 'domain-alignment'}."
        )
        raise ValueError


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

        match config.seq_kind:

            case "nucleotide":
                if not query_result.nucleotide_sequence:
                    logger.warning(f"No nucleotide sequence found for query '{query_result.query_id}'.")
                    continue
                fasta = Fasta(
                    header=query_result.query_id,
                    sequence=query_result.nucleotide_sequence if query_result.nucleotide_sequence else ""
                )

            case "protein":
                if not query_result.aminoacid_sequence:
                    logger.warning(f"No protein sequence found for query '{query_result.query_id}'.")
                    continue
                fasta = Fasta(
                    header=query_result.query_id,
                    sequence=query_result.aminoacid_sequence if query_result.aminoacid_sequence else ""
                )

            case "domain-alignment":
                if not query_result.domain_hits[0].domain_alignments[0].alignment_fragments[0].sequence:
                    logger.warning(f"No domain alignment found for query '{query_result.query_id}'.")
                    continue
                fasta = Fasta(
                    header=query_result.query_id,
                    sequence=query_result.domain_hits[0].domain_alignments[0].alignment_fragments[0].sequence if query_result.domain_hits[0].domain_alignments[0].alignment_fragments[0].sequence else ""
                )

            case _:
                logger.error(f"Sequence kind '{config.seq_kind}' not recognized.")
                sys.exit(1)


        try:
            print(fasta)
        except BrokenPipeError:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())
            sys.exit(1)


if __name__ == "__main__":
    run(sys.argv[1:])
