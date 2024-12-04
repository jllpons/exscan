#!/usr/bin/env python

r"""
Merge bedtools intersect output with serialized query results.

Usage:
    domtblop.py gffintersect [options] <serialized_domtblout> <bedtools_intersect_output>

Arguments:
    <serialized_domtblout> STR      Path to the serialized domtblout file.
                                    Use "-" to read from stdin. (default: "-")
    <bedtools_intersect_output>     The output of bedtools intersect with the -wa -wb flags.
                                    See below for an example.


Options:
    -h, --help                      Show this help message and exit
    -l, --loglevel STR              Set the logging level [default: INFO]

Examples:

    $ domtblop.py togff <domtblout> > domtblout.gff3 \
        && bedtools intersect -wa -wb -a domtblout.gff3 -b path/to/gff_db/* > intersect_output.txt \
        && domtblop.py gffintersect <serialized_domtblout> intersect_output.txt

    ...or in one line:
    $ domtblop.py gffintersect <serialized_domtblout> <(bedtools intersect -wa -wb -a <(domtblop.py togff <domtblout>) -b path/to/gff_db/* -filenames)

Notes:

    TODO: Explain how things are matched
"""

import argparse
from dataclasses import dataclass
import json
import logging
import os
import sys
from typing import (
        List,
        Tuple,
        )

from domtblop_parser import (
        CustomEncoder,
        GffFeature,
        HmmscanQueryResult,
        UnexpectedQueryIdFormat,
        )
from domtblop_utils import (
        setup_logger,
        read_input,
    )


@dataclass
class BedtoolsIntersectOutputRow:
    """
    Dataclass to represent the output of bedtools intersect.

    Attributes:
        a: The first feature in the intersect.
        b: The second feature in the intersect.
    """

    a: GffFeature
    filename_b: str
    b: GffFeature


    @classmethod
    def from_bedtools_intersect(cls, line: str) -> "BedtoolsIntersectOutputRow":
        """
        Create a BedtoolsIntersectOutputRow instance from a bedtools intersect line.

        Args:
            line (str): A bedtools intersect line.

        Returns:
            BedtoolsIntersectOutputRow: A BedtoolsIntersectOutputRow instance.
        """
        fields = line.strip().split("\t")
        gff_a = fields[0:9]
        filename_b = fields[9]
        gff_b = fields[10:]
        return cls(
            a=GffFeature.from_gff3("\t".join(gff_a)),
            filename_b=filename_b,
            b=GffFeature.from_gff3("\t".join(gff_b)),
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
        "bedtools_intersect_output",
        metavar="<bedtools_intersect_output>",
        type=str,
        nargs="?",
        )

    # Optional arguments
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
            f"File {config.serialized_domtblout} "
            "does not exist."
        )
        raise FileNotFoundError

    if not config.bedtools_intersect_output:
        logger.error("No bedtools intersect output provided.")
        raise ValueError

    if not os.path.exists(config.bedtools_intersect_output) and config.bedtools_intersect_output != "-":
        logger.error(
            f"File {config.bedtools_intersect_output} "
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
        + ", ".join(f"{k}={v}" for k,v in config.__dict__.items())
        )

    try:
        serialized_domtblout = read_input(config.serialized_domtblout)
        bedtools_intersect_file_handle = read_input(config.bedtools_intersect_output)
    except FileNotFoundError:
        sys.exit(1)


    # Iterate over the serialized domtblout file and parse each line as a HmmscanQueryResult
    query_results = []
    for line in serialized_domtblout:
        logger.debug(f"Parsing serialized query result: {line}")
        try:
            query_result = HmmscanQueryResult.from_json(json.loads(line))
        except UnexpectedQueryIdFormat as e:
            logger.error(e)
            sys.exit(1)
        logger.debug(f"Succesfully parsed query result: {query_result}")
        query_results.append(query_result)


    # Create a hashmap to quickly find the index of a query result by its query_id
    query_results_index_hashmap = dict()
    for i, query_result in enumerate(query_results):
        query_results_index_hashmap[query_result.query_id.replace('"', '')] = i


    # Iterate over the bedtools intersect output and add the intersecting features to the corresponding serialized query result
    for line in bedtools_intersect_file_handle:
        logger.debug(f"Processing line: {line}")

        try:
            intersect = BedtoolsIntersectOutputRow.from_bedtools_intersect(line)
        except Exception as e:
            logger.error(f"Error parsing line: {line}")
            logger.error(e)
            sys.exit(1)

        intersect_a_id = intersect.a.attributes.get("ID")
        intersect_a_id = intersect_a_id.replace('"', '') if intersect_a_id else None
        if not intersect_a_id:
            logger.error(f"Feature A does not have an ID in attributes: {intersect.a}")
            logger.error(f"Problematic line: {line}")
            sys.exit(1)
        elif intersect_a_id not in query_results_index_hashmap:
            logger.warning(f"ID {intersect_a_id} not found in query results.")
            logger.warning(f"Problematic line: {line}")
            continue

        i = query_results_index_hashmap[intersect_a_id]
        if query_results[i].gff_intersecting_features is None:
            query_results[i].gff_intersecting_features = []

        query_results[query_results_index_hashmap[intersect_a_id]].gff_intersecting_features.append(intersect.b)


    # Print the updated query results
    for query_result in query_results:
        try:
            print(json.dumps(query_result, cls=CustomEncoder))
        except BrokenPipeError:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())
            sys.exit(1)

