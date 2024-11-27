#!/usr/bin/env python

"""
Plot informative plots from the serialized domtblout file.

Usage:
    domtblop.py plot [options] <serialized_domtblout>


Arguments:
    <serialized_domtblout> STR      Path to the serialized domtblout file.
                                    Use "-" to read from stdin. (default: "-")

Options:
    -h, --help                      Show this help message and exit
    -l, --loglevel STR              Set the logging level [default: INFO]
"""

import argparse
import json
import logging
import os
import sys
from typing import List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from domtblop_parser import (
    HmmscanQueryResult,
    UnexpectedQueryIdFormat,
)
from domtblop_utils import setup_logger, read_input


sns.set_theme()


pt = 1./72.27 # 72.27 points per inch

two_columns_figure = 539.*pt

my_width = two_columns_figure

golden_ratio = (1 + 5**0.5) / 2

HORIZONTAL_FIGSIZE = (my_width, my_width / golden_ratio)
VERTICAL_FIGSIZE = (my_width*golden_ratio, my_width)


def mk_plots_dir() -> None:
    """
    Make a directory for the plots to be saved in.
    """
    cwd = os.getcwd()
    if not os.path.exists(f"{cwd}/plots"):
        os.mkdir(f"{cwd}/plots")


def mk_nhits_per_profile_barplot(
        query_results: List[HmmscanQueryResult],
        logger: logging.Logger,
        ) -> None:
    """
    Make a barplot of the number of hits per profile.
    Barplot should be horizontal and sorted by number of hits.
    """

    data = {"hmmProfile": []}
    for query_result in query_results:
        for domain_hit in query_result.domain_hits:
            data["hmmProfile"].append(domain_hit.name)

    df = pd.DataFrame(data)

    f, ax = plt.subplots(figsize=HORIZONTAL_FIGSIZE)

    sns.set_color_codes("muted")

    sns.countplot(
        y="hmmProfile",
        data=df,
        ax=ax,
        order=df["hmmProfile"].value_counts().index,
        )

    ax.bar_label(ax.containers[0], fontsize=8)


    plt.xlabel("Number of hits", size=9)
    plt.ylabel("HMM Profile", size=9)
    plt.yticks(size=8,)
    ax.set_title("Number of hits per HMM profile")

    plt.savefig("plots/nhits_per_profile_barplot.svg", format="svg")

    plt.clf()
    plt.close()


def mk_nhits_per_profile_and_per_chromosome_barplot(
        query_results: List[HmmscanQueryResult],
        logger: logging.Logger,
        ) -> None:
    """
    """

    data = {"hmmProfile": [], "chromosome": []}
    for query_result in query_results:
        for domain_hit in query_result.domain_hits:
            data["hmmProfile"].append(domain_hit.name)
            data["chromosome"].append(query_result.parent_sequence.sequence_id)

    df = pd.DataFrame(data)

    f, ax = plt.subplots(figsize=VERTICAL_FIGSIZE)

    sns.set_color_codes("muted")

    sns.countplot(
        data=df,
        x="chromosome",
        ax=ax,
        hue="hmmProfile",
        hue_order=df["hmmProfile"].value_counts().index,
        )

    ax.bar_label(ax.containers[0], fontsize=8)
    for container in ax.containers:
        ax.bar_label(container, fontsize=8)

    plt.xlabel("Chromosome", size=9)
    plt.ylabel("Number of hits", size=9)
    plt.xticks(size=8,)
    ax.set_title("Number of hits per HMM profile and per chromosome")

    plt.legend(loc="upper right", fontsize=8)

    plt.savefig("plots/nhits_per_profile_and_per_chromosome_barplot.svg", format="svg")



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

    query_results = []
    for line in file_handle:

        try:
            query_result = HmmscanQueryResult.from_json(json.loads(line))
        except UnexpectedQueryIdFormat as e:
            logger.error(e)
            sys.exit(1)

        query_results.append(query_result)


    mk_plots_dir()

    mk_nhits_per_profile_barplot(query_results, logger)

    mk_nhits_per_profile_and_per_chromosome_barplot(query_results, logger)
