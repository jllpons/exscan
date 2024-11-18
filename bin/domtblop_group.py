#!/usr/bin/env python

"""
Group query results from a serialized domtblout file by distance.

Usage:
    domtblop.py group [options] <serialized_domtblout>


Arguments:
    <serialized_domtblout> STR      Path to the serialized domtblout file.
                                    Use "-" to read from stdin. (default: "-")

    -d, --distance INT              Distance in base pairs to group sequences by.

Options:
    --both-strands                  Group sequences by distance on both strands. [default: False]
    -h, --help                      Show this help message and exit
    -l, --loglevel STR               Set the logging level [default: INFO]
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
        HmmscanQueryResult,
        UnexpectedQueryIdFormat,
        CustomEncoder,
        GroupedQueryResults,
        )
from domtblop_utils import (
        setup_logger,
        read_input,
    )


@dataclass
class QueryGroup:
    """
    Dataclass to hold a group of query results.
    """
    query_results: List[HmmscanQueryResult]


    def start(self) -> int:
        """
        Returns the start position of the first sequence in the group.

        Returns:
            int: Start position of the first sequence in the group
        """

        if not self.query_results:
            raise ValueError("Attempted to get start position of empty group.")

        return min([result.parent_sequence.start for result in self.query_results])


    def end(self) -> int:
        """
        Returns the end position of the last sequence in the group.

        Returns:
            int: End position of the last sequence in the group
        """

        if not self.query_results:
            raise ValueError("Attempted to get end position of empty group.")

        return max([result.parent_sequence.end for result in self.query_results])


    def group_id(self) -> str:
        """
        Returns the group id.

        Returns:
            str: Group id
        """

        return f"{self.start()}..{self.end()}"


    def n_hits_pos_strand(self) -> int:
        """
        """

        return len(
                [i for i in self.query_results if i.parent_sequence.strand == "+"]
                )


    def n_hits_neg_strand(self) -> int:
        """
        """

        return len(
                [i for i in self.query_results if i.parent_sequence.strand == "-"]
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
        "-d", "--distance",
        type=int,
        required=True,
    )

    # Optional arguments
    parser.add_argument(
        "--both-strands",
        action="store_true",
        default=False,
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


    if config.distance < 0:
        logger.error("Negative integers are not allowed for distance parameter")
        raise ValueError

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
        file_handle = read_input(config.serialized_domtblout)
    except FileNotFoundError:
        sys.exit(1)


    logger.debug("Reading and storing all Query Results into memory")
    query_results = []
    for line in file_handle:

        try:
            query_result = HmmscanQueryResult.from_json(json.loads(line))
        except UnexpectedQueryIdFormat as e:
            logger.error(e)
            sys.exit(1)

        query_results.append(query_result)
    logger.debug(f"Total number of parsed Query Results: {len(query_results)}")


    logger.debug("Grouping Query Results by distance")
    # If we sort all query results by start pos, we avoid multiple iterations
    query_results.sort(key=lambda x: x.parent_sequence.start)
    groups = []
    for query_result in query_results:
        logger.debug(f"Attempting to group Query Result with id: {query_result.query_id}")

        if not groups:
            groups.append(QueryGroup([query_result]))
            continue

        prev_group = groups[-1]
        given_distance = prev_group.end() + config.distance

        # Since query results have been sorted by start position, we can
        # safely assume current query result does not stat before the previous one
        if query_result.parent_sequence.start < given_distance:
            logger.debug(
                f"Query Result with id: {query_result.query_id} has start "
                + f"position {query_result.parent_sequence.start} which is less "
                + "than the sum of the previous group's end position "
                + f"({prev_group.end()}) and the distance ({config.distance})."
                + f" Given distance: {given_distance}. Current Query Result "
                + f"end position minus given distance: {query_result.parent_sequence.start - given_distance}"
                + f". Grouping with previous group: {prev_group.group_id()}"
                )
            groups[-1].query_results.append(query_result)
            continue

        logger.debug(
            f"Query Result with id: {query_result.query_id} has start "
            + f"position {query_result.parent_sequence.start} which is greater "
            + "than the sum of the previous group's end position "
            + f"({prev_group.end()}) and the distance ({config.distance})."
            + f" Given distance: {given_distance}. Current Query Results "
            + f"end position minus given distance: {query_result.parent_sequence.start - given_distance}"
            + f". Creating new group: {given_distance}..{query_result.parent_sequence.start}"
            )
        groups.append(QueryGroup([query_result]))
    logger.debug(f"Total number of groups: {len(groups)}")
    logger.debug([group.group_id() for group in groups])


    for group in groups:
        for query_result in group.query_results:

            query_result.group = GroupedQueryResults(
                    group_id=group.group_id(),
                    start=group.start(),
                    end=group.end(),
                    n_hits_pos_strand=group.n_hits_pos_strand(),
                    n_hits_neg_strand=group.n_hits_neg_strand(),
                    )

            try:
                print(json.dumps(query_result, cls=CustomEncoder))

            except BrokenPipeError:
                devnull = os.open(os.devnull, os.O_WRONLY)
                os.dup2(devnull, sys.stdout.fileno())
                sys.exit(1)

