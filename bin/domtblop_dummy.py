#!/usr/bin/env python

"""
domtblop.py dummy: Print a serialized dummy hit for testing purposes.

Usage: domtblop.py dummy

Options:
    -l --loglevel   Log level (default: INFO).
    -h --help       Show this screen.
    -v --version    Show version.
"""

import argparse
import json
import logging
import sys
from typing import Tuple, List


from domtblop_utils import (
        setup_logger,
    )
from domtblop_parser import (
        CustomEncoder,
        HmmscanQueryResult,
        HmmscanDomainHit,
        HmmscanDomainAlignment,
        DomainAlignmentFragment,
        GroupedQueryResults,
        UnexpectedQueryIdFormat
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
        default=False,
    )

    return parser


def setup_config(args: List) -> Tuple[argparse.Namespace, logging.Logger]:
    """
    Setup configuration for the script.

    Returns:
    argparse.Namespace: Parsed command-line arguments.
    """
    parser = setup_argparse()
    config = parser.parse_args(args)

    if config.help:
        print(__doc__)
        sys.exit(0)

    logger = setup_logger(config.loglevel)


    return config, logger


def build_dummy_hit() -> HmmscanQueryResult:
    """
    Build a dummy hit for testing purposes.

    Returns:
        Hit: Dummy

    Raises:
        UnexpectedHitIdFormat: If the hit ID format is unexpected.
    """

    dummy_query_result = HmmscanQueryResult(
            query_id="myFavoriteGene_frame=1_begin=1_end=100",
            sequence="MAGNIFICENTSEQUENCEAMAGNIFICENTSEQUENCEAMAGNIFICENTSEQUENCEAMAGNIFICENTSEQUENCEAMAGNIFICENTSEQUENCEA",
            domain_hits=[
                HmmscanDomainHit(
                    accession="PF00001.1",
                    name="PF00001",
                    description="Such a very interesting Domain",
                    full_sequence_evalue=1e-10,
                    full_sequence_score=100,
                    bias=0.9,
                    domain_alignments=[
                        HmmscanDomainAlignment(
                            independent_evalue=1e-07,
                            conditional_evalue=1e-05,
                            average_alignment_accuracy=0.9,
                            bias=0.5,
                            bit_score=50,
                            domain_number=1,
                            envelope_start=1,
                            envelope_end=30,
                            alignment_fragments=[
                                DomainAlignmentFragment(
                                    domain_start=1,
                                    domain_end=30,
                                    domain_strand="+",
                                    sequence_start=1,
                                    sequence_end=30,
                                    sequence_strand="+",
                                    ),
                                ],
                            ),
                        HmmscanDomainAlignment(
                            independent_evalue=1e-02,
                            conditional_evalue=0.01,
                            average_alignment_accuracy=0.9,
                            bias=0.5,
                            bit_score=50,
                            domain_number=2,
                            envelope_start=1,
                            envelope_end=30,
                            alignment_fragments=[
                                DomainAlignmentFragment(
                                    domain_start=1,
                                    domain_end=30,
                                    domain_strand="+",
                                    sequence_start=60,
                                    sequence_end=90,
                                    sequence_strand="+",
                                    ),
                                ],
                            ),

                        ]
                    ),
                HmmscanDomainHit(
                    accession="PF00002.1",
                    name="PF00002",
                    description="Not that interesting Domain",
                    full_sequence_evalue=1e-02,
                    full_sequence_score=17,
                    bias=2,
                    domain_alignments=[
                        HmmscanDomainAlignment(
                            independent_evalue=0.3,
                            conditional_evalue=0.01,
                            average_alignment_accuracy=0.9,
                            bias=0.5,
                            bit_score=3,
                            domain_number=1,
                            envelope_start=15,
                            envelope_end=25,
                            alignment_fragments=[
                                DomainAlignmentFragment(
                                    domain_start=8,
                                    domain_end=15,
                                    domain_strand="+",
                                    sequence_start=3,
                                    sequence_end=10,
                                    sequence_strand="+",
                                    ),
                                ],
                            ),
                        ],
                    ),
                ],
                group=GroupedQueryResults(
                    group_id="1..100",
                    start=1,
                    end=100,
                    n_hits_pos_strand=1,
                    n_hits_neg_strand=0
                )
            )

    dummy_query_result.compute_domain_aligment_fragment_absolut_positions()

    return dummy_query_result


def run(args: List) -> None:

    try:
        _, logger = setup_config(args)
    except (ValueError, FileNotFoundError):
        sys.exit(1)

    try:
        dummy_hit = build_dummy_hit()
    except UnexpectedQueryIdFormat:
        logger.error("Ironically, the dummy hit ID format is unexpected. Please fix this.")
        sys.exit(1)

    print(json.dumps(dummy_hit, cls=CustomEncoder))

