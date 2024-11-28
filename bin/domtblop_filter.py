#!/usr/bin/env python

"""
Apply filters to serialized domtblout query results.

Usage:
    domtblop.py filter [options] <serialized_domtblout>


Arguments:
    <serialized_domtblout> STR              Path to the serialized domtblout file.
                                            Use "-" to read from stdin. (default: "-")

Filtering options:

    Query-Sequence Level Filters:
    Those filter apply to the different hits a query sequence may have against different profile HMMs.

        --seq-evalue <value> FLOAT          Keep hits with full sequence e-value <= value
        --seq-score <value> FLOAT           Keep hits with full sequence score >= value
        --seq-bias <value> FLOAT            Keep hits with full sequence bias <= value

    Domain-Hit Level Filters:
    Those filter apply to the different hits a query sequence may have against the same profile HMM.

        --dom-ievalue <value> FLOAT         Keep hits with independent e-value <= value
        --dom-cevalue <value> FLOAT         Keep hits with conditional e-value <= value
        --dom-score <value> FLOAT           Keep hits with domain score >= value
        --dom-bias <value> FLOAT            Keep hits with domain bias <= value

    Domain-Alignment-Level Filters:
    Those filter apply to the different domain alignments a query sequence may have against the same profile HMM.

        --min-alignment-len <value> INT     Keep hits with domain alignment length >= value (sequence side)

Selection Options:

    --best-hit BOOL                         Only keep the best hit for each sequence or domain
                                            **based on filtering criteria**

Options:
    -h, --help                              Show this help message and exit
    -l, --loglevel STR                      Set the logging level [default: INFO]

Notes:

- **Comparasion Operators**:
    - For e-values, and bias options: values **EQUAL OR LESS THAN** the provided value are kept.
    - For scores: values **EQUAL OR GREATER THAN** the provided value are kept.

Examples:
    Filter by domain independent e-value, keeping only the best domain alignment:
      $ domtblop.py filter --dom-ievalue 1e-10 --best-hit-each-domain < cat query_results.json

    Filter the results with a full sequence E-value equal or smaller than 1e-05
      $ domtblop.py filter --seq-evalue 1e-05 < query_results.json

    For each query result, keep only the best hit based on the full sequence score
      $ domtblop.py filter --seq-score 0 --best-hit < query_results.json

    Only keep domain alignemnts if at least 50 residues are aligned
      $ domtblop.py filter --min-alignment-length 50 < query_results.json
"""

import argparse
from enum import Enum
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
        HmmscanDomainAlignment,
        HmmscanDomainHit,
        HmmscanQueryResult,
        UnexpectedQueryIdFormat,
        )
from domtblop_utils import (
        setup_logger,
        read_input,
    )


class QuerySequenceLevelFilter(Enum):
    EVALUE = "evalue"
    SCORE = "score"
    BIAS = "bias"


class DomainHitLevelFilter(Enum):
    IEVALUE = "ievalue"
    CEVALUE = "cevalue"
    SCORE = "score"
    BIAS = "bias"


class DomainAlignmentLevelFilter(Enum):
    ALIGNMENT_LENGTH = "alignment_length"


class QueryResultWithEmptyFields(Exception):
    pass


def filter_query_sequence_level(
        query_result: HmmscanQueryResult,
        value: int | float,
        filter_type: QuerySequenceLevelFilter,
        logger: logging.Logger,
        ) -> None:
    """
    Filter the query result at the query sequence level.

    Args:
        query_result (HmmscanQueryResult): Query result to filter.
        value (int | float): Value to filter against.
        filter_type (QuerySequenceLevelFilter): Type of filter to apply.
        logger (logging.Logger): Logger instance.

    Returns:
        None

    Raises:
        QueryResultWithEmptyFields: If the query result has no domain hits.

    Modifies:
        query_result: Removes domain hits that do not meet the filtering criteria.
    """
    logger.debug(
        f"Filtering query result {query_result.query_id} at the query sequence level"
        )

    if not query_result.domain_hits:
        raise QueryResultWithEmptyFields(
                f"Query result {query_result.query_id} has no domain hits"
                )

    match filter_type:

        case QuerySequenceLevelFilter.EVALUE:
            query_result.domain_hits = [
                dh for dh in query_result.domain_hits if dh.full_sequence_evalue <= float(value)
                ]

        case QuerySequenceLevelFilter.SCORE:
            query_result.domain_hits = [
                dh for dh in query_result.domain_hits if dh.full_sequence_score >= int(value)
                ]

        case QuerySequenceLevelFilter.BIAS:
            query_result.domain_hits = [
                dh for dh in query_result.domain_hits if dh.bias <= int(value)
                ]


def has_domain_hits(query_result: HmmscanQueryResult) -> bool:
    """
    Check if the query result has domain hits.

    Args:
        query_result (HmmscanQueryResult): Query result to check.

    Returns:
        bool: True if the query result has domain hits, False otherwise.
    """
    return any(query_result.domain_hits)


def keep_best_domain_hit(
        query_result: HmmscanQueryResult,
        filter_type: QuerySequenceLevelFilter,
        logger: logging.Logger,
        ) -> None:
    """
    Keep only the best domain hit based on the filtering criteria.

    Args:
        query_result (HmmscanQueryResult): Query result to filter.
        filter_type (QuerySequenceLevelFilter): Type of filter to apply.
        logger (logging.Logger): Logger instance.

    Returns:
        None

    Raises:
        QueryResultWithEmptyFields: If the query result has no domain hits.

    Modifies:
        query_result: Removes domain hits that do not meet the filtering criteria.
    """
    logger.debug(
        f"Keeping the best domain hit for query result {query_result.query_id}"
        )

    if not query_result.domain_hits:
        raise QueryResultWithEmptyFields(
                f"Query result {query_result.query_id} has no domain hits"
                )

    match filter_type:

        case QuerySequenceLevelFilter.EVALUE:
            best_hit = min(query_result.domain_hits, key=lambda x: x.full_sequence_evalue)
            query_result.domain_hits = [best_hit]

        case QuerySequenceLevelFilter.SCORE:
            best_hit = max(query_result.domain_hits, key=lambda x: x.full_sequence_score)
            query_result.domain_hits = [best_hit]

        case QuerySequenceLevelFilter.BIAS:
            best_hit = min(query_result.domain_hits, key=lambda x: x.bias)
            query_result.domain_hits = [best_hit]


def filter_domain_hit_level(
        query_result: HmmscanQueryResult,
        value: int | float,
        filter_type: DomainHitLevelFilter,
        logger: logging.Logger,
        ) -> None:
    """
    Filter the query result at the domain hit level.

    Args:
        query_result (HmmscanQueryResult): Query result to filter.
        value (int | float): Value to filter against.
        filter_type (DomainHitLevelFilter): Type of filter to apply.
        logger (logging.Logger): Logger instance.

    Returns:
        None

    Raises:
        QueryResultWithEmptyFields: If the query result has no domain hits.

    Modifies:
        query_result: Removes domain hits that do not meet the filtering criteria.
    """
    logger.debug(
        f"Filtering query result {query_result.query_id} at the domain hit level"
        )

    if not query_result.domain_hits:
        raise QueryResultWithEmptyFields(
                f"Query result {query_result.query_id} has no domain hits"
                )

    for domain_hit in query_result.domain_hits:

        if not domain_hit.domain_alignments:
            raise QueryResultWithEmptyFields(
                f"Query result {query_result.query_id} has no domain alignments"
                )

        match filter_type:

            case DomainHitLevelFilter.IEVALUE:
                domain_hit.domain_alignments = [
                    da for da in domain_hit.domain_alignments if da.independent_evalue <= float(value)
                    ]

            case DomainHitLevelFilter.CEVALUE:
                domain_hit.domain_alignments = [
                    da for da in domain_hit.domain_alignments if da.conditional_evalue <= float(value)
                    ]

            case DomainHitLevelFilter.SCORE:
                domain_hit.domain_alignments = [
                    da for da in domain_hit.domain_alignments if da.bit_score >= int(value)
                    ]

            case DomainHitLevelFilter.BIAS:
                domain_hit.domain_alignments = [
                    da for da in domain_hit.domain_alignments if da.bias <= int(value)
                    ]


def keep_best_domain_alignment(
        query_result: HmmscanQueryResult,
        filter_type: DomainHitLevelFilter,
        logger: logging.Logger,
        ) -> None:
    """
    Keep only the best domain alignment based on the filtering criteria.

    Args:
        query_result (HmmscanQueryResult): Query result to filter.
        filter_type (DomainHitLevelFilter): Type of filter to apply.
        logger (logging.Logger): Logger instance.

    Returns:
        None

    Raises:
        QueryResultWithEmptyFields: If the query result has no domain hits.

    Modifies:
        query_result: Removes domain hits that do not meet the filtering criteria.
    """

    logger.debug(
        f"Keeping the best domain alignment for query result {query_result.query_id}"
        )

    if not query_result.domain_hits:
        raise QueryResultWithEmptyFields(
                f"Query result {query_result.query_id} has no domain hits"
                )

    for domain_hit in query_result.domain_hits:

        if not domain_hit.domain_alignments:
            raise QueryResultWithEmptyFields(
                f"Query result {query_result.query_id} has no domain alignments"
                )

        match filter_type:

            case DomainHitLevelFilter.IEVALUE:
                best_alignment = min(domain_hit.domain_alignments, key=lambda x: x.independent_evalue)
                domain_hit.domain_alignments = [best_alignment]

            case DomainHitLevelFilter.CEVALUE:
                best_alignment = min(domain_hit.domain_alignments, key=lambda x: x.conditional_evalue)
                domain_hit.domain_alignments = [best_alignment]

            case DomainHitLevelFilter.SCORE:
                best_alignment = max(domain_hit.domain_alignments, key=lambda x: x.bit_score)
                domain_hit.domain_alignments = [best_alignment]

            case DomainHitLevelFilter.BIAS:
                best_alignment = min(domain_hit.domain_alignments, key=lambda x: x.bias)
                domain_hit.domain_alignments = [best_alignment]


def has_domain_alignments(domain_hit: HmmscanDomainHit) -> bool:
    """
    Check if the domain hit has domain alignments.

    Args:
        domain_hit (HmmscanDomainHit): Domain hit to check.

    Returns:
        bool: True if the domain hit has domain alignments, False otherwise.
    """

    return any(domain_hit.domain_alignments)


def filter_domain_alignment_level(
        query_result: HmmscanQueryResult,
        value: int,
        filter_type: DomainAlignmentLevelFilter,
        logger: logging.Logger,
        ) -> None:
    """
    Filter the domain hit at the domain alignment level.

    Args:
        query_result (HmmscanQueryResult): Query result to filter.
        value (int): Value to filter against.
        filter_type (DomainAlignmentLevelFilter): Type of filter to apply.
        logger (logging.Logger): Logger instance.

    Returns:
        None

    Raises:
        QueryResultWithEmptyFields: If the domain hit has no domain alignments.

    Modifies:
        domain_hit: Removes domain alignments that do not meet the filtering criteria.
    """

    logger.debug(
        f"Filtering query result {query_result.query_id} at the domain alignment level"
        )

    if not query_result.domain_hits:
        raise QueryResultWithEmptyFields(
                f"Query result {query_result.query_id} has no domain hits"
                )

    for domain_hit in query_result.domain_hits:

        if not domain_hit.domain_alignments:
            raise QueryResultWithEmptyFields(
                f"Query result {query_result.query_id} has no domain alignments"
                )


        for domain_alignment in domain_hit.domain_alignments:

            if not domain_alignment.alignment_fragments:
                raise QueryResultWithEmptyFields(
                    f"Query result {query_result.query_id} has no domain alignment fragments"
                    )

            match filter_type:

                case DomainAlignmentLevelFilter.ALIGNMENT_LENGTH:

                    if len(domain_alignment.alignment_fragments) == 1:
                        if domain_alignment.alignment_fragments[0].sequence_end - domain_alignment.alignment_fragments[0].sequence_start < value:
                            domain_hit.domain_alignments = []

                    else:
                        alignment_length = sum(fragment.sequence_end - fragment.sequence_start for fragment in domain_alignment.alignment_fragments)

                        if alignment_length < value:
                            domain_hit.domain_alignments = []


def has_alignment_fragments(domain_alignment: HmmscanDomainAlignment) -> bool:
    """
    Check if the domain alignment has alignment fragments.

    Args:
        domain_alignment (HmmscanDomainAlignment): Domain alignment to check.

    Returns:
        bool: True if the domain alignment has alignment fragments, False otherwise.
    """

    return any(domain_alignment.alignment_fragments)


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

    # Filtering options

    # Full sequence filters
    seq_filters = parser.add_mutually_exclusive_group()
    seq_filters.add_argument(
        "--seq-evalue",
        type=float,
        default=None,
        help="Keep hits with full sequence e-value < value"
    )
    seq_filters.add_argument(
        "--seq-score",
        type=float,
        default=None,
        help="Keep hits with full sequence score > value"
    )
    seq_filters.add_argument(
        "--seq-bias",
        type=float,
        default=None,
        help="Keep hits with full sequence bias < value"
    )

    # Domain filters
    dom_filters = parser.add_mutually_exclusive_group()
    dom_filters.add_argument(
        "--dom-ievalue",
        type=float,
        default=None,
        help="Keep hits with independent e-value < value"
    )
    dom_filters.add_argument(
        "--dom-cevalue",
        type=float,
        default=None,
        help="Keep hits with conditional e-value < value"
    )
    dom_filters.add_argument(
        "--dom-score",
        type=float,
        default=None,
        help="Keep hits with domain score > value"
    )
    dom_filters.add_argument(
        "--dom-bias",
        type=float,
        default=None,
        help="Keep hits with domain bias < value"
    )

    # Domain alignment filters
    parser.add_argument(
        "--min-alignment-length",
        type=int,
        default=None,
        help="Keep hits with domain alignment length > value (sequence side)"
    )

    # Optional arguments
    parser.add_argument(
        "--best-hit",
        action="store_true",
        default=False,
        help="Only keep the best hit for each sequence based on filtering criteria"
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
            f"Serialized domtblout file '{config.serialized_domtblout}' "
            "does not exist."
        )
        raise FileNotFoundError


    if not any([
        config.seq_evalue,
        config.seq_score,
        config.seq_bias,
        config.dom_ievalue,
        config.dom_cevalue,
        config.dom_score,
        config.dom_bias,
        config.min_alignment_length,
        ]):
        logger.error(
            "No filtering criteria provided. "
            "Please provide at least one filtering criterion."
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
        + ", ".join(f"{k}={v}" for k,v in config.__dict__.items())
        )

    try:
        file_handle = read_input(config.serialized_domtblout)
    except FileNotFoundError:
        sys.exit(1)

    n_hits_before_filter = 0
    n_hits_after_filter = 0
    n_domain_alignments_before_filter = 0
    n_domain_alignments_after_filter = 0
    for line in file_handle:

        try:
            query_result = HmmscanQueryResult.from_json(json.loads(line))
        except UnexpectedQueryIdFormat as e:
            logger.error(e)
            sys.exit(1)

        logger.debug(f"Attempting to apply filtering conditions to {query_result.query_id}")

        n_hits_before_filter += len(query_result.domain_hits)
        n_domain_alignments_before_filter += sum(len(domain_hit.domain_alignments) for domain_hit in query_result.domain_hits)

        if not any(query_result.domain_hits):
            logger.error(
                f"No domain hit where found for the query: {query_result.query_id}"
            )
            sys.exit(1)

        # Here we're applying the filters to whole sequence hit results
        if config.seq_evalue:
            filter_query_sequence_level(query_result, config.seq_evalue, QuerySequenceLevelFilter.EVALUE, logger)
            if config.best_hit:
                keep_best_domain_hit(query_result, QuerySequenceLevelFilter.EVALUE, logger)

        elif config.seq_score:
            filter_query_sequence_level(query_result, config.seq_score, QuerySequenceLevelFilter.SCORE, logger)
            if config.best_hit:
                keep_best_domain_hit(query_result, QuerySequenceLevelFilter.SCORE, logger)

        elif config.seq_bias:
            filter_query_sequence_level(query_result, config.seq_bias, QuerySequenceLevelFilter.BIAS, logger)
            if config.best_hit:
                keep_best_domain_hit(query_result, QuerySequenceLevelFilter.BIAS, logger)


        # Here we're applying the filters to domain hits
        if config.dom_ievalue:
            filter_domain_hit_level(query_result, config.dom_ievalue, DomainHitLevelFilter.IEVALUE, logger)
            if config.best_hit:
                keep_best_domain_alignment(query_result, DomainHitLevelFilter.IEVALUE, logger)

        elif config.dom_cevalue:
            filter_domain_hit_level(query_result, config.dom_cevalue, DomainHitLevelFilter.CEVALUE, logger)
            if config.best_hit:
                keep_best_domain_alignment(query_result, DomainHitLevelFilter.CEVALUE, logger)

        elif config.dom_score:
            filter_domain_hit_level(query_result, config.dom_score, DomainHitLevelFilter.SCORE, logger)
            if config.best_hit:
                keep_best_domain_alignment(query_result, DomainHitLevelFilter.SCORE, logger)

        elif config.dom_bias:
            filter_domain_hit_level(query_result, config.dom_bias, DomainHitLevelFilter.BIAS, logger)
            if config.best_hit:
                keep_best_domain_alignment(query_result, DomainHitLevelFilter.BIAS, logger)


        # Here we're applying the filters to domain alignments
        if config.min_alignment_length:
            filter_domain_alignment_level(query_result, config.min_alignment_length, DomainAlignmentLevelFilter.ALIGNMENT_LENGTH, logger)


        # NOTE: Remember the following structure:
        # query_result: [
        #           domain_hits: [
        #               domain_alignments: [
        #                   alignment_fragments: []
        #               ]
        #           ]
        #       ]

        # Has any hits after filtering?
        # query_result.domain_hits = [?]
        if not has_domain_hits(query_result):
            # Query result has no domain hits left, do not print it and skip to the next query result
            # rm(query_result.domain_hits = [])
            logger.warning(
                "After applying filters, there are no remaining domain hits for the query: "
                f"{query_result.query_id}"
                )
            continue

        # Has any domain alignments after filtering?
        # query_result.domain_hits = [domain_hit1.domain_alignments[?], ...]
        for domain_hit in query_result.domain_hits:
            if not has_domain_alignments(domain_hit):
                # Remove the domain hit if there are no domain alignments left
                # query_result.domain_hits = [rm(domain_hit1.domain_alignments[]), ...]
                query_result.domain_hits.remove(domain_hit)
                logger.warning(
                    "After applying filters, there are no remaining domain alignments for the query: "
                    f"{query_result.query_id}"
                    )
        # Check again if there are any domain hits left after possible removals
        # query_result.domain_hits = [?]
        if not any(query_result.domain_hits):
            # Query result has no domain hits left, do not print it and skip to the next query result
            # rm(query_result.domain_hits = [])
            logger.warning(
                "After applying filters, there are no remaining domain hits for the query: "
                f"{query_result.query_id}"
                )
            continue


        # Has any alignment fragments after filtering?
        # query_result.domain_hits = [domain_hit1.domain_alignments[domain_alignment1.alignment_fragments[?], ...], ...]
        for domain_hit in query_result.domain_hits:
            for domain_alignment in domain_hit.domain_alignments:
                if not has_alignment_fragments(domain_alignment):
                    # Remove the domain alignment if there are no alignment fragments left
                    # query_result.domain_hits = [domain_hit1.domain_alignments[rm(domain_alignment1.alignment_fragments[]), ...], ...]
                    domain_hit.domain_alignments.remove(domain_alignment)
                    logger.warning(
                        "After applying filters, there are no remaining domain alignment fragments for the query: "
                        f"{query_result.query_id}"
                        )
            if not has_domain_alignments(domain_hit):
                # Remove the domain hit if there are no domain alignments left
                # query_result.domain_hits = [rm(domain_hit1.domain_alignments[]), ...]
                query_result.domain_hits.remove(domain_hit)
                logger.warning(
                    "After applying filters, there are no remaining domain alignments for the query: "
                    f"{query_result.query_id}"
                    )
        if not any(query_result.domain_hits):
            # Query result has no domain hits left, do not print it and skip to the next query result
            # rm(query_result.domain_hits = [])
            logger.warning(
                "After applying filters, there are no remaining domain hits for the query: "
                f"{query_result.query_id}"
                )
            continue


        n_hits_after_filter += len(query_result.domain_hits)
        n_domain_alignments_after_filter += sum(len(domain_hit.domain_alignments) for domain_hit in query_result.domain_hits)
        try:

            print(json.dumps(query_result, cls=CustomEncoder))

        except BrokenPipeError:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())
            sys.exit(1)

    logger.info(
        f"Number of hits before filtering: {n_hits_before_filter}, "
        f"Number of hits after filtering: {n_hits_after_filter}, "
        f"Difference: {n_hits_before_filter - n_hits_after_filter}"
    )
    logger.info(
        f"Number of domain alignments before filtering: {n_domain_alignments_before_filter}, "
        f"Number of domain alignments after filtering: {n_domain_alignments_after_filter}, "
        f"Difference: {n_domain_alignments_before_filter - n_domain_alignments_after_filter}"
    )

