#!/usr/bin/env python

"""
Apply filters to serialized domtblout query results.

Usage:
    domtblop.py filter [options] <serialized_domtblout>


Arguments:
    <serialized_domtblout> STR      Path to the serialized domtblout file.
                                    Use "-" to read from stdin. (default: "-")

Filtering options:

    Sequence-Level Filters:
    Those filter apply to the different hits a query sequence may have against different profile HMMs.

        --seq-evalue <value> FLOAT     Keep hits with full sequence e-value < value
        --seq-score <value> FLOAT      Keep hits with full sequence score > value
        --seq-bias <value> FLOAT       Keep hits with full sequence bias < value

    Domain-Level Filters:
    Those filter apply to the different hits a query sequence may have against the same profile HMM.

        --dom-ievalue <value> FLOAT    Keep hits with independent e-value < value
        --dom-cevalue <value> FLOAT    Keep hits with conditional e-value < value
        --dom-score <value> FLOAT      Keep hits with domain score > value
        --dom-bias <value> FLOAT       Keep hits with domain bias < value

Selection Options:

    --best-hit BOOL                    Only keep the best hit for each sequence or domain
                                       **based on filtering criteria**

Options:
    -h, --help                         Show this help message and exit
    -l, --loglevel STR                 Set the logging level [default: INFO]

Notes:

- **Comparasion Operators**:
    - For e-values, and bias options: values **LESS THAN** the provided value are kept.
    - For scores: values **GREATER THAN** the provided value are kept.

Examples:
    Filter by domain independent e-value, keeping only the best domain alignment:
      $ domtblop.py filter --dom-ievalue 1e-10 --best-hit-each-domain < cat query_results.json

    Filter the results with a full sequence E-value smaller than 1e-05
      $ domtblop.py filter --seq-evalue 1e-05 < query_results.json

    For each query result, keep only the best hit based on the full sequence score
      $ domtblop.py filter --seq-score 0 --best-hit < query_results.json
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
        HmmscanQueryResult,
        UnexpectedQueryIdFormat,
        CustomEncoder
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


    if not any([
        config.seq_evalue,
        config.seq_score,
        config.seq_bias,
        config.dom_ievalue,
        config.dom_cevalue,
        config.dom_score,
        config.dom_bias,
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

    for line in file_handle:

        try:
            query_result = HmmscanQueryResult.from_json(json.loads(line))
        except UnexpectedQueryIdFormat as e:
            logger.error(e)
            sys.exit(1)

        logger.debug(f"Attempting to apply filtering conditions to {query_result.query_id}")

        if not any(query_result.domain_hits):
            # When applying filters to a newly parsed hmmscan domtblout file,
            # we should not expect any query sequences without domain hits.
            # But if different filters are chained, it is possible that
            # some query sequences have no domain hits left after applying
            # previous filters.
            logger.warning(
                f"No domain hit where found for the query: {query_result.query_id}"
            )
            continue

        # Here we're applying the filters to whole sequence hit results
        if config.seq_evalue:
            query_result.filter_hits(config.seq_evalue, "evalue")
            if config.best_hit:
                query_result.keep_best_domain_hit("evalue")


        if config.seq_score:
            query_result.filter_hits(config.seq_score, "score")
            if config.best_hit:
                query_result.keep_best_domain_hit("score")


        if config.seq_bias:
            query_result.filter_hits(config.seq_bias, "bias")
            if config.best_hit:
                query_result.keep_best_domain_hit("bias")


        # Here we're applying the filters to domain hits
        for domain_hit in query_result.domain_hits:
            if config.dom_ievalue:
                domain_hit.filter_alignments(config.dom_ievalue, "ievalue")
                if config.best_hit:
                    domain_hit.keep_best_domain_alignment("ievalue")

            if config.dom_cevalue:
                domain_hit.filter_alignments(config.dom_cevalue, "cevalue")
                if config.best_hit:
                    domain_hit.keep_best_domain_alignment("cevalue")

            if config.dom_score:
                domain_hit.filter_alignments(config.dom_score, "score")
                if config.best_hit:
                    domain_hit.keep_best_domain_alignment("score")

            if config.dom_bias:
                domain_hit.filter_alignments(config.dom_bias, "bias")
                if config.best_hit:
                    domain_hit.keep_best_domain_alignment("bias")

        try:
            print(json.dumps(query_result, cls=CustomEncoder))

        except BrokenPipeError:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())
            sys.exit(1)



