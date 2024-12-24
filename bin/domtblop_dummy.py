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
        UnexpectedQueryIdFormat,
        CustomEncoder,
        GffFeature,
        GroupedQueryResults,
        DomainAlignmentFragment,
        HmmscanQueryResult,
        HmmscanDomainHit,
        HmmscanDomainAlignment,
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
            query_id="myFavoriteGene_frame=2_begin=10_end=200",
            aminoacid_sequence="MSEQNFFICENTSEQUENCEMSEQNFFICENTSEQUENCEMSEQNFFICENTSEQUENCEMSEQNFFICENTSEQUENCE",
            nucleotide_sequence="ATGAGTGAACAGNAATTTCATNNNTNNTTNNNNNNNNNNNNAATTTCTEEESOMEDATA... (truncated for brevity)",
            domain_hits=[
                HmmscanDomainHit(
                    accession="PF00010.5",
                    name="Interestingdomain",
                    description="A highly interesting domain for testin",
                    full_sequence_evalue=1e-12,
                    full_sequence_score=150,
                    bias=0.1,
                    domain_alignments=[
                        HmmscanDomainAlignment(
                            independent_evalue=1e-15,
                            conditional_evalue=1e-12,
                            average_alignment_accuracy=0.95,
                            bias=0.05,
                            bit_score=80,
                            domain_number=1,
                            envelope_start=10,
                            envelope_end=60,
                            alignment_fragments=[
                                DomainAlignmentFragment(
                                    domain_start=1,
                                    domain_end=25,
                                    domain_strand="+",
                                    sequence_start=19,
                                    sequence_end=94,
                                    sequence_strand="+",
                                    ),
                                DomainAlignmentFragment(
                                    domain_start=26,
                                    domain_end=30,
                                    domain_strand="+",
                                    sequence_start=36,
                                    sequence_end=42,
                                    sequence_strand="+",
                                    ),
                                ],
                            ),
                        HmmscanDomainAlignment(
                            independent_evalue=0.0001,
                            conditional_evalue=0.0001,
                            average_alignment_accuracy=0.88,
                            bias=0.5,
                            bit_score=45,
                            domain_number=2,
                            envelope_start=50,
                            envelope_end=90,
                            alignment_fragments=[
                                DomainAlignmentFragment(
                                    domain_start=1,
                                    domain_end=15,
                                    domain_strand="+",
                                    sequence_start=50,
                                    sequence_end=64,
                                    sequence_strand="+",
                                    ),
                                ],
                            ),

                        ]
                    ),
                HmmscanDomainHit(
                    accession="PF00011.3",
                    name="BorderlineDomain",
                    description="A domain with borderline values",
                    full_sequence_evalue=1e-05,
                    full_sequence_score=20,
                    bias=2,
                    domain_alignments=[
                        HmmscanDomainAlignment(
                            independent_evalue=0.01,
                            conditional_evalue=0.05,
                            average_alignment_accuracy=0.9,
                            bias=1.0,
                            bit_score=10,
                            domain_number=1,
                            envelope_start=70,
                            envelope_end=100,
                            alignment_fragments=[
                                DomainAlignmentFragment(
                                    domain_start=5,
                                    domain_end=14,
                                    domain_strand="+",
                                    sequence_start=75,
                                    sequence_end=84,
                                    sequence_strand="+",
                                    ),
                                ],
                            ),
                        ],
                    ),
                HmmscanDomainHit(
                    accession="PF00012.1",
                    name="LowEvalueDomain",
                    description="A domain that likely won't pass strict filters",
                    full_sequence_evalue=0.1,
                    full_sequence_score=5,
                    bias=5.0,
                    domain_alignments=[
                        HmmscanDomainAlignment(
                            independent_evalue=0.5,
                            conditional_evalue=0.5,
                            average_alignment_accuracy=0.7,
                            bias=2.5,
                            bit_score=2,
                            domain_number=1,
                            envelope_start=90,
                            envelope_end=110,
                            alignment_fragments=[
                                DomainAlignmentFragment(
                                    domain_start=1,
                                    domain_end=5,
                                    domain_strand="+",
                                    sequence_start=95,
                                    sequence_end=99,
                                    sequence_strand="+",
                                    ),
                                ],
                            ),
                        ],
                    ),
                HmmscanDomainHit(
                    accession="PF00014.3",
                    name="ShortAlignmentDomain",
                    description="A domain alignment that is very short",
                    full_sequence_evalue=1e-08,
                    full_sequence_score=40,
                    bias=0.8,
                    domain_alignments=[
                        HmmscanDomainAlignment(
                            independent_evalue=1e-06,
                            conditional_evalue=1e-06,
                            average_alignment_accuracy=0.93,
                            bias=0.3,
                            bit_score=25,
                            domain_number=1,
                            envelope_start=30,
                            envelope_end=33,
                            alignment_fragments=[
                                DomainAlignmentFragment(
                                    domain_start=2,
                                    domain_end=3,
                                    domain_strand="+",
                                    sequence_start=31,
                                    sequence_end=32,
                                    sequence_strand="+",
                                    ),
                                ],
                            ),
                        ],
                    ),
                ],
                gff_intersecting_features=[
                    GffFeature(
                        seqid="ExampleGene",
                        source="trustedSource",
                        type_="gene",
                        start=1,
                        end=1000,
                        score=".",
                        strand="+",
                        phase=".",
                        attributes={"ID": "ExampleGene", "Name": "ExampleGene"}
                        ),
                    GffFeature(
                        seqid="ExampleGene",
                        source="trustedSource",
                        type_="exon",
                        start=20,
                        end=50,
                        score=".",
                        strand="+",
                        phase=".",
                        attributes={"ID": "ExampleGene.exon1", "Parent": "ExampleGene"}
                        ),
                    GffFeature(
                        seqid="ExampleGene",
                        source="trustedSource",
                        type_="CDS",
                        start=60,
                        end=63,
                        score=".",
                        strand="+",
                        phase=".",
                        attributes={"ID": "ExampleGene.CDS.1", "Parent": "ExampleGene"}
                        ),
                ],
                group=GroupedQueryResults(
                    group_id="10..200",
                    start=10,
                    end=200,
                    n_hits_pos_strand=3,
                    n_hits_neg_strand=0
                )
            )

    dummy_query_result.compute_domain_aligment_fragment_absolut_positions()
    dummy_query_result.mk_domain_alignment_fragment_sequences()

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

