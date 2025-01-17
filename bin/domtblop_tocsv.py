#!/usr/bin/env python

"""
Store the information of the serialized query results in a CSV file.

Usage:
    domtblop.py tocsv [options] <serialized_domtblout>


Arguments:
    <serialized_domtblout> STR      Path to the serialized domtblout file.
                                    Use "-" to read from stdin. (default: "-")

Options:
    --kind STR                      The kind of output to generate.
                                    Either sequence-centric or domain alignment-centric.
                                    Choices: ["sequence", "domain"] [default: "sequence"]

    --separator STR                 The separator to use in the CSV file. [default: ","]
    -h, --help                      Show this help message and exit
    -l, --loglevel STR              Set the logging level [default: INFO]

Examples:
    $ domtblop.py tocsv domtblout.json
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
    Set,
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

@dataclass
class SequenceCentricTableRecord:
    """
    _domain_hits: Dict[str, int]  Binary representation of the domain hits.
                  {
                      "domain1": 1,
                      "domain2": 0,
                      "domain3": 1,
                      "domain4": 0,
                  }
    """

    query_id: str
    #sequence_type: str
    # parent_sequence_id actually goes here in the csv output
    _domain_hits: Dict[str, int]
    # group_id: str
    # gff_features: List[GffFeature]
    # protein_sequence: str
    # source_sequence: str

    parent_sequence_id: str | None = None
    #group_id: str | None = None
    #gff_features: List[GffFeature] | None = None
    protein_sequence: str | None = None
    source_sequence: str | None = None

    def __post_init__(self):
        self.protein_length = str(len(self.protein_sequence)) if self.protein_sequence else ""

    def to_csv(self, separator: str = ",") -> str:

        domain_hits = [
                str(observation)
                for _domain, observation in self._domain_hits.items()
        ]


        return separator.join(
            [
                self.query_id,
                #self.sequence_type,
                self.parent_sequence_id if self.parent_sequence_id else "",
                separator.join(domain_hits),
                #self.group_id if self.group_id else "",
                #"::".join([f'"{str(gff.to_json())}"' for gff in self.gff_features]) if self.gff_features else "",
                self.protein_length,
                self.protein_sequence if self.protein_sequence else "",
                self.source_sequence if self.source_sequence else "",
            ]
        )

    def csv_headers(self, separator: str = ",") -> str:


        domain_hits = [
            domain
            for domain, _observation in self._domain_hits.items()
        ]

        return separator.join(
            [
                "query_id",
                #"sequence_type",
                "parent_sequence_id",
                separator.join(domain_hits),
                #"group_id",
                #"gff_features",
                "protein_length",
                "protein_sequence",
                "source_sequence",
            ]
        )

    @classmethod
    def from_query_result(
        cls,
        query_result: HmmscanQueryResult,
        domain_set: Set[str],
        ) -> "SequenceCentricTableRecord":

        domain_hits = {domain: 0 for domain in domain_set}
        for domain_hit in query_result.domain_hits:
            if domain_hit.name in domain_set:
                domain_hits[domain_hit.name] = 1

        return cls(
            query_id=query_result.query_id,
            #sequence_type=query_result.sequence_type,
            parent_sequence_id=query_result.parent_sequence.sequence_id if query_result.parent_sequence else None,
            _domain_hits=domain_hits,
            #group_id=query_result.group.id if query_result.group else None,
            #gff_features=query_result.gff_features,
            protein_sequence=query_result.protein_sequence,
            source_sequence=query_result.source_sequence,
        )


@dataclass
class DomainAlignmentCentricTableRecord:

    query_id: str
    domain_accession: str
    domain_name: str
    independent_evalue: float
    conditional_evalue: float
    average_alignment_accuracy: float
    bias: float
    bitscore: float
    domain_start: int
    domain_end: int
    sequence_start: int
    sequence_end: int
    protein_sequence: str | None = None
    gff_features: List[GffFeature] | None = None

    def __post_init__(self):
        self.alignment_length = str(self.sequence_end - self.sequence_start + 1)

    def to_csv(self, separator: str = ",") -> str:

        return separator.join(
            [
                self.query_id,
                self.domain_accession,
                self.domain_name,
                str(self.independent_evalue),
                str(self.conditional_evalue),
                str(self.average_alignment_accuracy),
                str(self.bias),
                str(self.bitscore),
                str(self.domain_start),
                str(self.domain_end),
                str(self.sequence_start),
                str(self.sequence_end),
                self.alignment_length,
                self.protein_sequence if self.protein_sequence else "",
                ";".join([gff.feature_id for gff in self.gff_features]) if self.gff_features else "",
            ]
        )

    def csv_headers(self, separator: str = ",") -> str:
        return separator.join(
            [
                "query_id",
                "domain_accession",
                "domain_name",
                "independent_evalue",
                "conditional_evalue",
                "average_alignment_accuracy",
                "bias",
                "bitscore",
                "domain_start",
                "domain_end",
                "sequence_start",
                "sequence_end",
                "alignment_length",
                "protein_sequence",
                "gff_features",
            ]
        )

    @classmethod
    def from_domain_hit(
        cls,
        query_result: HmmscanQueryResult
        ) -> List["DomainAlignmentCentricTableRecord"]:

        records = []

        for domain_hit in query_result.domain_hits:
            for domain_alignment in domain_hit.domain_alignments:
                for domain_alignment_fragment in domain_alignment.alignment_fragments:
                    records.append(
                        cls(
                            query_id=query_result.query_id,
                            domain_accession=domain_hit.accession,
                            domain_name=domain_hit.name,
                            independent_evalue=domain_alignment.independent_evalue,
                            conditional_evalue=domain_alignment.conditional_evalue,
                            average_alignment_accuracy=domain_alignment.average_alignment_accuracy,
                            bias=domain_alignment.bias,
                            bitscore=domain_alignment.bit_score,
                            domain_start=domain_alignment_fragment.domain_start,
                            domain_end=domain_alignment_fragment.domain_end,
                            sequence_start=domain_alignment_fragment.sequence_start,
                            sequence_end=domain_alignment_fragment.sequence_end,
                            protein_sequence=query_result.protein_sequence,
                            gff_features=query_result.gff_features,
                        )
                    )

        return records


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
        "--kind",
        type=str,
        default="sequence",
        choices=["sequence", "domain"],
    )
    parser.add_argument(
        "--separator",
        type=str,
        default=",",
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

    query_results = []
    for line in file_handle:
        try:
            query_result = HmmscanQueryResult.from_json(json.loads(line))
        except UnexpectedQueryIdFormat as e:
            logger.error(e)
            sys.exit(1)

        query_results.append(query_result)

    domain_set = set()
    for query_result in query_results:
        for domain_hit in query_result.domain_hits:
            domain_set.add(domain_hit.name)

    table_records = []
    match config.kind:
        case "domain":
            table_records.extend(
                record
                for query_result in query_results
                for record in DomainAlignmentCentricTableRecord.from_domain_hit(query_result)
            )

        case "sequence":
            table_records = [
                SequenceCentricTableRecord.from_query_result(query_result, domain_set)
                for query_result in query_results
            ]

        case _:
            logger.error(f"Invalid kind of output: {config.kind}")

    for i, record in enumerate(table_records):
        try:
            if i == 0:
                print(record.csv_headers(config.separator))

            print(record.to_csv(config.separator))

        except BrokenPipeError:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())
            sys.exit(1)
