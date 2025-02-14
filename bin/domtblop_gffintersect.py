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
    Dict,
    List,
    Tuple,
)

from domtblop_parser import (
    CustomEncoder,
    GffFeature,
    HmmscanQueryResult,
    Intersection,
    UnexpectedQueryIdFormat,
)
from domtblop_utils import (
    setup_logger,
    read_input,
)


@dataclass
class GffintersectParams:
    bedtools_intersect_output: str
    name: str = "gffintersect"

    def to_json(self) -> Dict[str, str]:
        return {
            "name": self.name,
            "bedtools_intersect_output": self.bedtools_intersect_output,
        }


@dataclass
class BedtoolsIntersectOutputRow:
    """
    Dataclass to represent the output of bedtools intersect.

    Attributes:
        a: The first feature in the intersect.
        b: The second feature in the intersect.
    """

    a: GffFeature
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
        gff_b = fields[9:]
        return cls(
            a=GffFeature.from_gff3("\t".join(gff_a)),
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

    if (
        not os.path.exists(config.serialized_domtblout)
        and config.serialized_domtblout != "-"
    ):
        logger.error(f"File {config.serialized_domtblout} " "does not exist.")
        raise FileNotFoundError

    if not config.bedtools_intersect_output:
        logger.error("No bedtools intersect output provided.")
        raise ValueError

    if (
        not os.path.exists(config.bedtools_intersect_output)
        and config.bedtools_intersect_output != "-"
    ):
        logger.error(f"File {config.bedtools_intersect_output} " "does not exist.")
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

    # Given a list of query results, create a hashmap where the key is the query_id
    # and the value is the index of the query result in the list.
    # This allows us to quickly find the index of a query result in the list by its query_id
    query_results_index_hashmap = dict()
    for i, query_result in enumerate(query_results):
        query_results_index_hashmap[query_result.query_id.replace('"', "")] = i

    # Iterate over the bedtools intersect output and add the intersecting features to the corresponding serialized query result
    # 'intersect_b_feature_hashmap' is used to keep track of the unique feature in the intersect_b column
    # Keep in mind that a feature in the intersect_b column can be repeated multiple times,
    # where each time will intersect with a different hmm domain alignment.
    # We want to asign a unique feature_id to each unique feature in the intersect_b column
    # so we can reference it in the serialized query result.
    intersect_b_feature_hashmap = dict()
    n_unique_features = 1
    for line in bedtools_intersect_file_handle:
        logger.debug(f"Processing line: {line}")

        try:
            intersect = BedtoolsIntersectOutputRow.from_bedtools_intersect(line)
        except Exception as e:
            logger.error(f"Error parsing line: {line}")
            logger.error(e)
            sys.exit(1)

        if intersect.b.to_gff3() not in intersect_b_feature_hashmap:
            intersect_b_feature_hashmap[intersect.b.to_gff3()] = n_unique_features
            n_unique_features += 1

        intersect.b.feature_id = (
            f"GF_{str(intersect_b_feature_hashmap.get(intersect.b.to_gff3()))}"
        )

        intersect_a_id = intersect.a.attributes.get("parentID", None)
        intersect_a_id = intersect_a_id.replace('"', "") if intersect_a_id else None
        if not intersect_a_id:
            logger.error(
                f"Feature A does not have an parentID in attributes: {intersect.a}"
            )
            logger.error(f"Problematic line: {line}")
            sys.exit(1)
        elif intersect_a_id not in query_results_index_hashmap:
            logger.warning(f"ID {intersect_a_id} not found in query results.")
            logger.warning(f"Problematic line: {line}")
            continue

        i = query_results_index_hashmap[intersect_a_id]
        query_result = query_results[i]
        if query_result.gff_features is None:
            query_result.gff_features = []
        if query_result.intersections is None:
            query_result.intersections = []

        if intersect.b.feature_id not in [
            f.feature_id for f in query_result.gff_features
        ]:
            query_result.gff_features.append(intersect.b)

        query_result.intersections.append(
            Intersection(
                domain_alignment_id=intersect.a.feature_id,
                feature_id=intersect.b.feature_id,
            )
        )

    # Print the updated query results
    for query_result in query_results:
        try:
            query_result.metadata.parameters_used.append(
                GffintersectParams(
                    bedtools_intersect_output=config.bedtools_intersect_output,
                )
            )
            query_result.update_modified_at()

            print(json.dumps(query_result, cls=CustomEncoder))
        except BrokenPipeError:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())
            sys.exit(1)
