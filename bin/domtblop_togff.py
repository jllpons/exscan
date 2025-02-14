#!/usr/bin/env python

"""
Write the profile HMM alignments from each serialized query result to stdout in GFF format.

Usage:
    domtblop.py togff [options] <serialized_domtblout>


Arguments:
    <serialized_domtblout> STR      Path to the serialized domtblout file.
                                    Use "-" to read from stdin. (default: "-")

Options:
    --intersecting-features         Output the intersecting of each query result instead
    -h, --help                      Show this help message and exit
    -l, --loglevel STR              Set the logging level [default: INFO]

Examples:
    $ domtblop.py togff domtblout.json > domtblout.gff3

    the gff features look like this:
seq1    domtblop.py     hmmscan_hit_alignment    1       100     1e-07     +       .       ID="query_id";DomainName="PF00001"
"""

import argparse
import json
import logging
import os
import sys
from typing import (
    List,
    Tuple,
)

from domtblop_parser import (
    GffFeature,
    HmmscanQueryResult,
    UnexpectedQueryIdFormat,
)
from domtblop_utils import (
    setup_logger,
    read_input,
)


def gff3_from_query_result(query_result: HmmscanQueryResult) -> List[GffFeature]:
    """
    Convert a HmmscanQueryResult object to a list of Gff3Feature objects.
    Remember that each query result:
        - Can have multiple domain hits
        - Each domain hit can have multiple domain alignments
        - Each domain alignment can have multiple alignment fragments

    Each alignment fragment will be represented as a Gff3Feature object.

    Args:
        query_result (HmmscanQueryResult): A HmmscanQueryResult object.

    Returns:
        List[Gff3Feature]: A list of Gff3Feature objects.
    """
    if not query_result.parent_sequence:
        raise ValueError(
            f"No parent sequence found in query result: {query_result.query_id}"
        )
    features = []

    seqid = query_result.parent_sequence.sequence_id
    source = "domtblop.py"
    strand = query_result.parent_sequence.strand
    phase = abs(query_result.parent_sequence.frame)

    for domain_hit in query_result.domain_hits:
        for domain_alignment in domain_hit.domain_alignments:
            for domain_alignment_fragment in domain_alignment.alignment_fragments:
                type_ = "hmmscan_hit_alignment"
                # NOTE: hmmscan tblout seems to be 1-based, so I'll keep it that way
                start = (
                    domain_alignment_fragment.sequence_start_in_parent
                    if domain_alignment_fragment.sequence_start_in_parent
                    else -1
                )
                end = (
                    domain_alignment_fragment.sequence_end_in_parent
                    if domain_alignment_fragment.sequence_end_in_parent
                    else -1
                )
                score = domain_alignment.independent_evalue
                attributes = {
                    "featureID": domain_alignment_fragment.domain_alignment_id,
                    "parentID": query_result.query_id,
                    "DomainName": domain_hit.name,
                }

                if any([start < 0, end < 0]):
                    raise ValueError(
                        f"No value was found for start or end in alignment fragment: {domain_alignment_fragment}"
                    )

                features.append(
                    GffFeature(
                        feature_id=domain_alignment_fragment.domain_alignment_id,
                        seqid=seqid,
                        source=source,
                        type_=type_,
                        start=start,
                        end=end,
                        score=score,
                        strand=strand,
                        phase=phase,
                        attributes=attributes,
                    )
                )

    return features


def gff3_from_query_result_intersecting(
    query_result: HmmscanQueryResult,
) -> List[GffFeature]:
    """ """
    features = []

    if not query_result.gff_features:
        return features

    for intersecting_gff in query_result.gff_features:
        features.append(
            GffFeature(
                feature_id=intersecting_gff.feature_id,
                seqid=intersecting_gff.seqid,
                source=intersecting_gff.source,
                type_=intersecting_gff.type_,
                start=intersecting_gff.start,
                end=intersecting_gff.end,
                score=intersecting_gff.score,
                strand=intersecting_gff.strand,
                phase=intersecting_gff.phase,
                attributes=intersecting_gff.attributes,
            )
        )

    return features


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

    # Optional arguments
    parser.add_argument(
        "--intersecting-features",
        action="store_true",
        default=False,
        help="Output the intersecting of each query result instead",
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

        # TODO: Handle groups?
        try:
            if config.intersecting_features:
                features = gff3_from_query_result_intersecting(query_result)

            else:
                features = gff3_from_query_result(query_result)

        except ValueError as e:
            logger.error(f"Error processing query result: {e}")
            sys.exit(1)

        for feature in features:
            try:
                print(feature.to_gff3())

            except BrokenPipeError:
                devnull = os.open(os.devnull, os.O_WRONLY)
                os.dup2(devnull, sys.stdout.fileno())
                sys.exit(1)
