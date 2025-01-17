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

    Interecting GFF Features Level Filters:
    Those filter apply to the different GFF features that intersect with the domain alignments of a query sequence.

        --keep-only-intersect               Keep only the GFF features that intersect with the domain alignments of a query sequence


Selection Options:

    --best-hit BOOL                         Only keep the best hit for each sequence or domain **based on filtering criteria**

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
from dataclasses import dataclass
from enum import Enum
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
    HmmscanQueryResult,
    UnexpectedQueryIdFormat,
)
from domtblop_utils import (
    setup_logger,
    read_input,
)


@dataclass
class FilterParams:
    _filter_enum: Enum
    value: int | float | None
    keep_best: bool
    name: str = "filter"

    @property
    def filter_type(self) -> str:
        match self._filter_enum:
            case QuerySequenceLevelFilter.EVALUE:
                return "query_sequence_evalue"

            case QuerySequenceLevelFilter.SCORE:
                return "query_sequence_score"

            case QuerySequenceLevelFilter.BIAS:
                return "query_sequence_bias"

            case DomainHitLevelFilter.IEVALUE:
                return "domain_hit_ievalue"

            case DomainHitLevelFilter.CEVALUE:
                return "domain_hit_cevalue"

            case DomainHitLevelFilter.SCORE:
                return "domain_hit_score"

            case DomainHitLevelFilter.BIAS:
                return "domain_hit_bias"

            case DomainAlignmentLevelFilter.ALIGNMENT_LENGTH:
                return "domain_alignment_length"

            case IntersectingGFFFeaturesLevelFilter.KEEP_ONLY_INTERSECT:
                return "keep_only_intersect"

            case _:
                raise ValueError(f"Invalid filter type: {self._filter_enum}")

    def to_json(self) -> Dict:
        return {
            "name": self.name,
            "filter_type": self.filter_type,
            "value": self.value,
            "keep_best": self.keep_best,
        }


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


class IntersectingGFFFeaturesLevelFilter(Enum):  # not used yet
    KEEP_ONLY_INTERSECT = "keep_only_intersect"


class QueryResultWithEmptyFields(Exception):
    pass


def validate_query_result(query_result: HmmscanQueryResult) -> bool:
    """
    Validate that query_result has the minimal structure required.
    Return True if valid, False if invalid.

    Args:
        query_result (HmmscanQueryResult): Query result to validate.

    Returns:
        bool: True if the query result has the minimal structure required, False otherwise.
    """
    if not query_result.domain_hits:
        return False
    for hit in query_result.domain_hits:
        if not hit.domain_alignments:
            return False
        for alignment in hit.domain_alignments:
            if not alignment.alignment_fragments:
                return False

    return True


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
                dh
                for dh in query_result.domain_hits
                if dh.full_sequence_evalue <= float(value)
            ]

        case QuerySequenceLevelFilter.SCORE:
            query_result.domain_hits = [
                dh
                for dh in query_result.domain_hits
                if dh.full_sequence_score >= int(value)
            ]

        case QuerySequenceLevelFilter.BIAS:
            query_result.domain_hits = [
                dh for dh in query_result.domain_hits if dh.bias <= float(value)
            ]


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
            best_hit = min(
                query_result.domain_hits, key=lambda x: x.full_sequence_evalue
            )
            query_result.domain_hits = [best_hit]

        case QuerySequenceLevelFilter.SCORE:
            best_hit = max(
                query_result.domain_hits, key=lambda x: x.full_sequence_score
            )
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
                    da
                    for da in domain_hit.domain_alignments
                    if da.independent_evalue <= float(value)
                ]

            case DomainHitLevelFilter.CEVALUE:
                domain_hit.domain_alignments = [
                    da
                    for da in domain_hit.domain_alignments
                    if da.conditional_evalue <= float(value)
                ]

            case DomainHitLevelFilter.SCORE:
                domain_hit.domain_alignments = [
                    da
                    for da in domain_hit.domain_alignments
                    if da.bit_score >= int(value)
                ]

            case DomainHitLevelFilter.BIAS:
                domain_hit.domain_alignments = [
                    da for da in domain_hit.domain_alignments if da.bias <= float(value)
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
                best_alignment = min(
                    domain_hit.domain_alignments, key=lambda x: x.independent_evalue
                )
                domain_hit.domain_alignments = [best_alignment]

            case DomainHitLevelFilter.CEVALUE:
                best_alignment = min(
                    domain_hit.domain_alignments, key=lambda x: x.conditional_evalue
                )
                domain_hit.domain_alignments = [best_alignment]

            case DomainHitLevelFilter.SCORE:
                best_alignment = max(
                    domain_hit.domain_alignments, key=lambda x: x.bit_score
                )
                domain_hit.domain_alignments = [best_alignment]

            case DomainHitLevelFilter.BIAS:
                best_alignment = min(domain_hit.domain_alignments, key=lambda x: x.bias)
                domain_hit.domain_alignments = [best_alignment]


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
                        if (
                            domain_alignment.alignment_fragments[0].sequence_end
                            - domain_alignment.alignment_fragments[0].sequence_start
                            < value
                        ):
                            domain_hit.domain_alignments = []

                    else:
                        alignment_length = sum(
                            fragment.sequence_end - fragment.sequence_start
                            for fragment in domain_alignment.alignment_fragments
                        )

                        if alignment_length < value:
                            domain_hit.domain_alignments = []


def cleanup_empty_structures(
    query_result: HmmscanQueryResult,
) -> HmmscanQueryResult | None:
    """
    Cleanup empty structures in the query result.

    Args:
        query_result (HmmscanQueryResult): Query result to cleanup.

    Returns:
        HmmscanQueryResult | None: Query result with empty structures removed or None if no data remains.
    """

    # Check if query result has domain hits. If not, function returns None
    if not query_result.has_domain_hits():
        return None

    # Iterate over hits to only retain only those that have alignments
    query_result.domain_hits = [
        dh for dh in query_result.domain_hits if dh.has_domain_alignments()
    ]

    # We may have removed domain hits in the previous operation, re-checking
    if not query_result.has_domain_hits():
        return None

    # Now we are left with domain hits that we know have domain alignments
    for dh in query_result.domain_hits:
        # So we retain only those domain alignments that contain alignment fragments
        dh.domain_alignments = [
            da for da in dh.domain_alignments if da.has_domain_alignment_fragments()
        ]

    # Retain only the domain hits that contain domain alginments
    query_result.domain_hits = [
        dh for dh in query_result.domain_hits if dh.has_domain_alignments()
    ]

    if not query_result.has_domain_hits():
        return None

    return query_result


def apply_filters(
    query_result: HmmscanQueryResult, config: argparse.Namespace, logger: logging.Logger
) -> HmmscanQueryResult | None:
    """
    Apply all relevant filters to the query_result based on config.
    Return the filtered query_result (possibly modified) or None if no data remains.

    Args:
        query_result (HmmscanQueryResult): Query result to filter.
        config (argparse.Namespace): Configuration options.
        logger (logging.Logger): Logger instance.

    Returns:
        HmmscanQueryResult | None: Filtered query result or None if no data remains.
    """
    # Apply sequence-level filters
    if config.seq_evalue is not None:
        filter_query_sequence_level(
            query_result, config.seq_evalue, QuerySequenceLevelFilter.EVALUE, logger
        )
        query_result = cleanup_empty_structures(query_result)
        if query_result is None:
            return None
        if config.best_hit:
            keep_best_domain_hit(query_result, QuerySequenceLevelFilter.EVALUE, logger)

        query_result.metadata.parameters_used.append(
            FilterParams(
                QuerySequenceLevelFilter.EVALUE,
                config.seq_evalue,
                config.best_hit,
            )
        )

    elif config.seq_score is not None:
        filter_query_sequence_level(
            query_result, config.seq_score, QuerySequenceLevelFilter.SCORE, logger
        )
        query_result = cleanup_empty_structures(query_result)
        if query_result is None:
            return None
        if config.best_hit:
            keep_best_domain_hit(query_result, QuerySequenceLevelFilter.SCORE, logger)

        query_result.metadata.parameters_used.append(
            FilterParams(
                QuerySequenceLevelFilter.SCORE,
                config.seq_score,
                config.best_hit,
            )
        )

    elif config.seq_bias is not None:
        filter_query_sequence_level(
            query_result, config.seq_bias, QuerySequenceLevelFilter.BIAS, logger
        )
        query_result = cleanup_empty_structures(query_result)
        if query_result is None:
            return None
        if config.best_hit:
            keep_best_domain_hit(query_result, QuerySequenceLevelFilter.BIAS, logger)

        query_result.metadata.parameters_used.append(
            FilterParams(
                QuerySequenceLevelFilter.BIAS,
                config.seq_bias,
                config.best_hit,
            )
        )

    # Domain-level filters
    if config.dom_ievalue is not None:
        filter_domain_hit_level(
            query_result, config.dom_ievalue, DomainHitLevelFilter.IEVALUE, logger
        )
        query_result = cleanup_empty_structures(query_result)
        if query_result is None:
            return None
        if config.best_hit:
            keep_best_domain_alignment(
                query_result, DomainHitLevelFilter.IEVALUE, logger
            )

        query_result.metadata.parameters_used.append(
            FilterParams(
                DomainHitLevelFilter.IEVALUE,
                config.dom_ievalue,
                config.best_hit,
            )
        )

    elif config.dom_cevalue is not None:
        filter_domain_hit_level(
            query_result, config.dom_cevalue, DomainHitLevelFilter.CEVALUE, logger
        )
        query_result = cleanup_empty_structures(query_result)
        if query_result is None:
            return None
        if config.best_hit:
            keep_best_domain_alignment(
                query_result, DomainHitLevelFilter.CEVALUE, logger
            )

        query_result.metadata.parameters_used.append(
            FilterParams(
                DomainHitLevelFilter.CEVALUE,
                config.dom_cevalue,
                config.best_hit,
            )
        )

    elif config.dom_score is not None:
        filter_domain_hit_level(
            query_result, config.dom_score, DomainHitLevelFilter.SCORE, logger
        )
        query_result = cleanup_empty_structures(query_result)
        if query_result is None:
            return None
        if config.best_hit:
            keep_best_domain_alignment(query_result, DomainHitLevelFilter.SCORE, logger)

        query_result.metadata.parameters_used.append(
            FilterParams(
                DomainHitLevelFilter.SCORE,
                config.dom_score,
                config.best_hit,
            )
        )

    elif config.dom_bias is not None:
        filter_domain_hit_level(
            query_result, config.dom_bias, DomainHitLevelFilter.BIAS, logger
        )
        query_result = cleanup_empty_structures(query_result)
        if query_result is None:
            return None
        if config.best_hit:
            keep_best_domain_alignment(query_result, DomainHitLevelFilter.BIAS, logger)

        query_result.metadata.parameters_used.append(
            FilterParams(
                DomainHitLevelFilter.BIAS,
                config.dom_bias,
                config.best_hit,
            )
        )

    # Here we're applying the filters to domain alignments:
    if config.min_alignment_length:
        filter_domain_alignment_level(
            query_result,
            config.min_alignment_length,
            DomainAlignmentLevelFilter.ALIGNMENT_LENGTH,
            logger,
        )
        query_result = cleanup_empty_structures(query_result)
        if query_result is None:
            return None

        query_result.metadata.parameters_used.append(
            FilterParams(
                DomainAlignmentLevelFilter.ALIGNMENT_LENGTH,
                config.min_alignment_length,
                False,
            )
        )

    # Cleanup empty structures after filtering
    query_result = cleanup_empty_structures(query_result)

    # If requested, keep only entries with intersecting GFF features
    if config.keep_only_intersect:
        query_result.metadata.parameters_used.append(
            FilterParams(
                IntersectingGFFFeaturesLevelFilter.KEEP_ONLY_INTERSECT,
                None,
                False,
            )
        )

        if not query_result.has_intersecting_gff_features():
            return None

    # If no hits remain after filtering and cleanup
    if not query_result.has_domain_hits():
        return None

    if not any(dh.has_domain_alignments() for dh in query_result.domain_hits):
        return None

    if not any(
        da.has_domain_alignment_fragments()
        for dh in query_result.domain_hits
        for da in dh.domain_alignments
    ):
        return None

    query_result.update_modified_at()

    return query_result


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
        help="Keep hits with full sequence e-value < value",
    )
    seq_filters.add_argument(
        "--seq-score",
        type=float,
        default=None,
        help="Keep hits with full sequence score > value",
    )
    seq_filters.add_argument(
        "--seq-bias",
        type=float,
        default=None,
        help="Keep hits with full sequence bias < value",
    )

    # Domain filters
    dom_filters = parser.add_mutually_exclusive_group()
    dom_filters.add_argument(
        "--dom-ievalue",
        type=float,
        default=None,
        help="Keep hits with independent e-value < value",
    )
    dom_filters.add_argument(
        "--dom-cevalue",
        type=float,
        default=None,
        help="Keep hits with conditional e-value < value",
    )
    dom_filters.add_argument(
        "--dom-score",
        type=float,
        default=None,
        help="Keep hits with domain score > value",
    )
    dom_filters.add_argument(
        "--dom-bias",
        type=float,
        default=None,
        help="Keep hits with domain bias < value",
    )

    # Domain alignment filters
    parser.add_argument(
        "--min-alignment-length",
        type=int,
        default=None,
        help="Keep hits with domain alignment length > value (sequence side)",
    )

    # Intersecting GFF features filters
    parser.add_argument(
        "--keep-only-intersect",
        action="store_true",
        default=False,
        help="Keep only the GFF features that intersect with the domain alignments of a query sequence",
    )

    # Optional arguments
    parser.add_argument(
        "--best-hit",
        action="store_true",
        default=False,
        help="Only keep the best hit for each sequence based on filtering criteria",
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
        logger.error(
            f"Serialized domtblout file '{config.serialized_domtblout}' "
            "does not exist."
        )
        raise FileNotFoundError

    if not any(
        [
            config.seq_evalue,
            config.seq_score,
            config.seq_bias,
            config.dom_ievalue,
            config.dom_cevalue,
            config.dom_score,
            config.dom_bias,
            config.min_alignment_length,
            config.keep_only_intersect,
        ]
    ):
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
        + ", ".join(f"{k}={v}" for k, v in config.__dict__.items())
    )

    try:
        file_handle = read_input(config.serialized_domtblout)
    except FileNotFoundError:
        sys.exit(1)

    n_qresults_before_filter = 0
    n_qresults_after_filter = 0
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

        logger.debug(
            f"Attempting to apply filtering conditions to {query_result.query_id}"
        )

        if not validate_query_result(query_result):
            logger.error(f"Invalid query result {query_result.query_id}")
            sys.exit(1)

        # Before filtering
        original_hits = len(query_result.domain_hits)
        original_alignments = sum(
            len(dh.domain_alignments) for dh in query_result.domain_hits
        )
        try:
            id_ = query_result.query_id
            query_result = apply_filters(query_result, config, logger)
        except QueryResultWithEmptyFields as e:
            logger.error(e)
            sys.exit(1)

        # If we are still here, we have a valid query result that we counted BEFORE filtering
        n_qresults_before_filter += 1
        n_hits_before_filter += original_hits
        n_domain_alignments_before_filter += original_alignments

        if query_result is None:
            logger.warning(f"Query result {id_} does not satisfy filtering criteria")
            continue

        n_qresults_after_filter += 1
        n_hits_after_filter += len(query_result.domain_hits)
        n_domain_alignments_after_filter += sum(
            len(dh.domain_alignments) for dh in query_result.domain_hits
        )

        try:
            print(json.dumps(query_result, cls=CustomEncoder))

        except BrokenPipeError:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())
            sys.exit(1)

    # Compute differences
    hits_removed = n_hits_before_filter - n_hits_after_filter
    alignments_removed = (
        n_domain_alignments_before_filter - n_domain_alignments_after_filter
    )

    # Handle division by zero (in case no queries processed)
    qresults_percentage = (
        (n_qresults_after_filter / n_qresults_before_filter * 100)
        if n_qresults_before_filter > 0
        else 0
    )
    hits_percentage = (
        (hits_removed / n_hits_before_filter * 100) if n_hits_before_filter > 0 else 0
    )
    alignments_percentage = (
        (alignments_removed / n_domain_alignments_before_filter * 100)
        if n_domain_alignments_before_filter > 0
        else 0
    )

    logger.info(
        "\nFiltering Summary:\n"
        "------------------\n"
        f"Query Results: {n_qresults_before_filter} before, {n_qresults_after_filter} after "
        f"({qresults_percentage:.2f}% remain)\n"
        f"Domain Hits: {n_hits_before_filter} before, {n_hits_after_filter} after, {hits_removed} removed "
        f"({hits_percentage:.2f}% removed)\n"
        f"Domain Alignments: {n_domain_alignments_before_filter} before, {n_domain_alignments_after_filter} after, "
        f"{alignments_removed} removed ({alignments_percentage:.2f}% removed)\n"
    )
