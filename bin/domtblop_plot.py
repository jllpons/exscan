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
import numpy as np
import pandas as pd
import seaborn as sns
from upsetplot import from_memberships, UpSet

from domtblop_parser import (
    HmmscanQueryResult,
    UnexpectedQueryIdFormat,
)
from domtblop_utils import setup_logger, read_input


# Set the default style for the plots
sns.set_theme()

# Set up the figure sizes
pt = 1.0 / 72.27  # 72.27 points per inch
two_columns_figure = 539.0 * pt
golden_ratio = (1 + 5**0.5) / 2

my_width = two_columns_figure

HORIZONTAL_FIGSIZE = (my_width, my_width / golden_ratio)
VERTICAL_FIGSIZE = (my_width * golden_ratio, my_width)
XXLARGE_HORIZ_FIGSIZE = (my_width * 5, my_width * 5 / golden_ratio)
XXLARGE_VERT_FIGSIZE = (my_width * 5 * golden_ratio, my_width * 5)


def mk_plots_dir() -> None:
    """
    Make a directory for the plots to be saved in.
    """
    cwd = os.getcwd()
    if not os.path.exists(f"{cwd}/plots"):
        os.mkdir(f"{cwd}/plots")


def mk_to_png_script() -> None:

    script = r"""
#!/usr/bin/env bash

# Convert all SVG files in the current directory to PNG
ls plots | grep ".svg" | sed 's/.svg//g' | xargs -I{}  resvg -w 3200 -h 2600 --dpi 300 plots/{}.svg plots/{}.png
    """

    with open(f"{os.getcwd()}/to_png.sh", "w") as f:
        f.write(script)


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
            data["hmmProfile"].append("Total")

    df = pd.DataFrame(data)

    f, ax = plt.subplots(figsize=HORIZONTAL_FIGSIZE)

    sns.set_color_codes("muted")

    sns.countplot(
        y="hmmProfile",
        data=df,
        ax=ax,
        order=df["hmmProfile"].value_counts().index,
    )

    try:
        ax.bar_label(ax.containers[0], fontsize=8)
    except IndexError:
        pass

    plt.xlabel("Number of hits", size=9)
    plt.ylabel("HMM Profile", size=9)
    plt.yticks(
        size=8,
    )
    ax.set_title("Number of hits per HMM profile")

    plt.savefig("plots/nhits_per_profile_barplot.svg", format="svg")

    plt.clf()
    plt.close()


def mk_nuniqhits_per_profile_sequence_level_barplot(
    query_results: List[HmmscanQueryResult],
    logger: logging.Logger,
    ) -> None:

    data = {"hmmProfile": []}
    for query_result in query_results:
        profile_hits = set()
        for domain_hit in query_result.domain_hits:
            profile_hits.add(domain_hit.name)
        data["hmmProfile"].extend(profile_hits)
        data["hmmProfile"].extend(["Total" for _ in range(len(profile_hits))])

    df = pd.DataFrame(data)

    f, ax = plt.subplots(figsize=HORIZONTAL_FIGSIZE)

    sns.set_color_codes("muted")

    sns.countplot(
        y="hmmProfile",
        data=df,
        ax=ax,
        order=df["hmmProfile"].value_counts().index,
        )

    try:
        ax.bar_label(ax.containers[0], fontsize=8)
    except IndexError:
        pass

    plt.xlabel("Number of hits", size=9)
    plt.ylabel("HMM Profile", size=9)
    plt.yticks(
        size=8,
    )
    ax.set_title("Number of unique hits per HMM profile for each sequence")

    plt.savefig("plots/nuniqhits_per_profile_sequence_level_barplot.svg", format="svg")

    plt.clf()
    plt.close()


def mk_nhits_per_profile_and_per_chromosome_barplot(
    query_results: List[HmmscanQueryResult],
    logger: logging.Logger,
) -> None:
    """ """

    data = {"hmmProfile": [], "chromosome": []}
    for query_result in query_results:
        for domain_hit in query_result.domain_hits:
            data["hmmProfile"].append(domain_hit.name)
            data["chromosome"].append(query_result.parent_sequence.sequence_id)

    df = pd.DataFrame(data)
    # Sort the df by chromosome (taking into account the chromosome number)
    df["chromosome_number"] = df["chromosome"].str.extract(r"(\d+)")
    df = df.sort_values(by=["chromosome_number", "chromosome", "hmmProfile"])

    f, ax = plt.subplots(figsize=XXLARGE_HORIZ_FIGSIZE)

    sns.set_color_codes("muted")

    sns.countplot(
        data=df,
        x="chromosome",
        ax=ax,
        hue="hmmProfile",
        hue_order=df["hmmProfile"].value_counts().index,
    )

    ax.bar_label(ax.containers[0], fontsize=9)
    for container in ax.containers:
        ax.bar_label(container, fontsize=9)

    plt.xlabel("Chromosome", size=30)
    plt.ylabel("Number of hits", size=30)
    plt.xticks(size=20, rotation=45)
    plt.yticks(
        size=20,
    )
    ax.set_title("Number of hits per HMM profile and per chromosome", size=40)

    plt.legend(loc="upper right", fontsize=20)

    plt.savefig("plots/nhits_per_profile_and_per_chromosome_barplot.svg", format="svg")
    plt.clf()
    plt.close()


def mk_violinplot_evalue_distribution_per_profile(
    query_results: List[HmmscanQueryResult],
    logger: logging.Logger,
) -> None:
    data = {"query_id": [], "hmmProfile": [], "evalue": []}
    for query_result in query_results:
        for domain_hit in query_result.domain_hits:
            data["query_id"].append(query_result.query_id)
            data["hmmProfile"].append(domain_hit.name)
            data["evalue"].append(domain_hit.full_sequence_evalue)

    df = pd.DataFrame(data)

    df = df.sort_values(by=["hmmProfile", "evalue"])

    f, ax = plt.subplots(figsize=HORIZONTAL_FIGSIZE)

    sns.violinplot(
        y="evalue",
        x="hmmProfile",
        data=df,
        ax=ax,
        scale="width",
        hue="hmmProfile",
    )

    # ax.set_yscale("log")
    # ax.set_ylim(min(df["evalue"]), 1e-3)
    plt.xlabel("HMM Profile", size=9)
    plt.ylabel("E-value", size=9)
    plt.yticks(
        size=8,
    )
    plt.xticks(size=6)
    ax.set_title("E-value distribution per HMM profile")

    plt.savefig("plots/violinplot_evalue_distribution_per_profile.svg", format="svg")
    plt.clf()
    plt.close()


def mk_violinplot_aligment_length_distribution_per_profile(
    query_results: List[HmmscanQueryResult],
    logger: logging.Logger,
) -> None:
    data = {"query_id": [], "hmmProfile": [], "alignment_length": []}
    for query_result in query_results:
        for domain_hit in query_result.domain_hits:
            data["query_id"].append(query_result.query_id)
            data["hmmProfile"].append(domain_hit.name)

            best_domain_alignment = max(
                domain_hit.domain_alignments, key=lambda x: x.independent_evalue
            )

            length = 0
            for domain_alignment_fragment in best_domain_alignment.alignment_fragments:
                length += (
                    domain_alignment_fragment.sequence_end
                    - domain_alignment_fragment.sequence_start
                )

            data["alignment_length"].append(length)

    df = pd.DataFrame(data)

    df.sort_values(by=["hmmProfile", "alignment_length"])

    f, ax = plt.subplots(figsize=HORIZONTAL_FIGSIZE)

    sns.violinplot(
        y="alignment_length",
        x="hmmProfile",
        data=df,
        ax=ax,
        scale="width",
        hue="hmmProfile",
    )

    plt.xlabel("HMM Profile", size=9)
    plt.ylabel("No. Residues", size=9)
    plt.yticks(
        size=8,
    )
    plt.xticks(size=6)
    ax.set_title(
        "No. of residues from translated sequence aligned against HMM profiles"
    )

    plt.savefig(
        "plots/violinplot_alignment_length_distribution_per_profile.svg", format="svg"
    )
    plt.clf()
    plt.close()


def mk_violin_bitscore_distribution_per_profile(
    query_results: List[HmmscanQueryResult],
    logger: logging.Logger,
) -> None:
    """ """

    data = {"query_id": [], "hmmProfile": [], "bitscore": []}
    for query_result in query_results:
        for domain_hit in query_result.domain_hits:
            data["query_id"].append(query_result.query_id)
            data["hmmProfile"].append(domain_hit.name)
            data["bitscore"].append(domain_hit.full_sequence_score)

    df = pd.DataFrame(data)

    df = df.sort_values(by=["hmmProfile", "bitscore"])

    f, ax = plt.subplots(figsize=HORIZONTAL_FIGSIZE)

    sns.violinplot(
        y="bitscore",
        x="hmmProfile",
        data=df,
        ax=ax,
        scale="width",
        hue="hmmProfile",
    )

    plt.xlabel("HMM Profile", size=9)
    plt.ylabel("Bitscore", size=9)

    plt.yticks(
        size=8,
    )
    plt.xticks(size=6)

    ax.set_title("Bitscore distribution per HMM profile Hits")

    plt.savefig("plots/violinplot_bitscore_distribution_per_profile.svg", format="svg")
    plt.clf()


def mk_upset_hits_per_sequence(
    query_results: List[HmmscanQueryResult],
    logger: logging.Logger,
) -> None:
    """ """

    query_id_memberships = {}
    for query_result in query_results:
        query_id_memberships[query_result.query_id] = ""

    for query_result in query_results:
        for domain_hit in query_result.domain_hits:
            if domain_hit.name not in query_id_memberships[query_result.query_id].split(
                ","
            ):
                query_id_memberships[query_result.query_id] += f"{domain_hit.name},"

    # remove trailing comma
    for query_id in query_id_memberships:
        query_id_memberships[query_id] = query_id_memberships[query_id][:-1]

    to_lists = [i.split(",") for i in query_id_memberships.values()]

    df = from_memberships(to_lists)

    f, ax = plt.subplots(figsize=XXLARGE_HORIZ_FIGSIZE)
    upset = UpSet(
        df,
        subset_size="count",
        sort_by="cardinality",
        show_counts=True,
        show_percentages=True,
        element_size=70,
    )

    upset.plot()

    plt.title("Profile HMM Hits per Query Sequence", size=25)

    plt.savefig("plots/upset_hits_per_qsequence.svg", format="svg")
    plt.clf()
    plt.close()


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

    mk_to_png_script()

    mk_nuniqhits_per_profile_sequence_level_barplot(query_results, logger)

    #mk_nhits_per_profile_barplot(query_results, logger)

    #mk_nhits_per_profile_and_per_chromosome_barplot(query_results, logger)

    mk_upset_hits_per_sequence(query_results, logger)

    #mk_violinplot_evalue_distribution_per_profile(query_results, logger)

    #mk_violinplot_aligment_length_distribution_per_profile(query_results, logger)

    #mk_violin_bitscore_distribution_per_profile(query_results, logger)
