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
                                    Choices: {"sequence", "domain", "gff", "map", "group"} [default: "sequence"]
                                      - "sequence": One row per sequence.
                                      - "domain": One row per domain alignment.
                                      - "gff": One row per each GFF feature that intersects with a domain alignment.
                                      - "map": Map intersecting GFF features to domain alignments.
                                      - "group": One row per group of domain alignments.

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
    # sequence_type: str
    # parent_sequence_id actually goes here in the csv output
    _domain_hits: Dict[str, int]
    # group_id: str
    # gff_features: List[GffFeature]
    # protein_sequence: str
    # source_sequence: str

    parent_sequence_id: str | None = None
    # group_id: str | None = None
    # gff_features: List[GffFeature] | None = None
    protein_sequence: str | None = None
    source_sequence: str | None = None

    def __post_init__(self):
        self.protein_length = (
            str(len(self.protein_sequence)) if self.protein_sequence else ""
        )

    def to_csv(self, separator: str = ",") -> str:
        domain_hits = [
            str(observation) for _domain, observation in self._domain_hits.items()
        ]

        return separator.join(
            [
                self.query_id,
                # self.sequence_type,
                self.parent_sequence_id if self.parent_sequence_id else "",
                separator.join(domain_hits),
                # self.group_id if self.group_id else "",
                # "::".join([f'"{str(gff.to_json())}"' for gff in self.gff_features]) if self.gff_features else "",
                self.protein_length,
                self.protein_sequence if self.protein_sequence else "",
                self.source_sequence if self.source_sequence else "",
            ]
        )

    def csv_headers(self, separator: str = ",") -> str:
        domain_hits = [domain for domain, _observation in self._domain_hits.items()]

        return separator.join(
            [
                "query_id",
                # "sequence_type",
                "parent_sequence_id",
                separator.join(domain_hits),
                # "group_id",
                # "gff_features",
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
            # sequence_type=query_result.sequence_type,
            parent_sequence_id=query_result.parent_sequence.sequence_id
            if query_result.parent_sequence
            else None,
            _domain_hits=domain_hits,
            # group_id=query_result.group.id if query_result.group else None,
            # gff_features=query_result.gff_features,
            protein_sequence=query_result.protein_sequence,
            source_sequence=query_result.source_sequence,
        )


@dataclass
class DomainAlignmentCentricTableRecord:
    domain_alignment_id: str
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
                self.domain_alignment_id,
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
                ";".join([gff.feature_id for gff in self.gff_features])
                if self.gff_features
                else "",
            ]
        )

    def csv_headers(self, separator: str = ",") -> str:
        return separator.join(
            [
                "domain_alignment_id",
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
        query_result: HmmscanQueryResult,
    ) -> List["DomainAlignmentCentricTableRecord"]:
        records = []

        for domain_hit in query_result.domain_hits:
            for domain_alignment in domain_hit.domain_alignments:
                for domain_alignment_fragment in domain_alignment.alignment_fragments:
                    records.append(
                        cls(
                            domain_alignment_id=domain_alignment_fragment.domain_alignment_id,
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


@dataclass
class IntersectingGffFeatureTableRecord:
    feature_id: str
    seqid: str
    source: str
    type_: str
    start: int
    end: int
    score: str | float
    strand: str
    phase: str | int
    attributes: Dict[str, str]

    def to_csv(self, separator: str = ",") -> str:
        return separator.join(
            [
                self.feature_id,
                self.seqid,
                self.source,
                self.type_,
                str(self.start),
                str(self.end),
                str(self.score),
                self.strand,
                str(self.phase),
                ";".join([f"{k}={v}" for k, v in self.attributes.items()]),
            ]
        )

    def csv_headers(self, separator: str = ",") -> str:
        return separator.join(
            [
                "feature_id",
                "seqid",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
            ]
        )

    @classmethod
    def from_query_result(
        cls,
        query_result: HmmscanQueryResult,
    ) -> List["IntersectingGffFeatureTableRecord"]:
        if not query_result.gff_features:
            return []

        records = []

        for gff_feature in query_result.gff_features:
            records.append(
                cls(
                    feature_id=gff_feature.feature_id,
                    seqid=gff_feature.seqid,
                    source=gff_feature.source,
                    type_=gff_feature.type_,
                    start=gff_feature.start,
                    end=gff_feature.end,
                    score=gff_feature.score,
                    strand=gff_feature.strand,
                    phase=gff_feature.phase,
                    attributes=gff_feature.attributes,
                )
            )

        return records


@dataclass
class IntersectionMapTableRecord:
    domain_alignment_id: str
    feature_id: str

    def to_csv(self, separator: str = ",") -> str:
        return separator.join(
            [
                self.domain_alignment_id,
                self.feature_id,
            ]
        )

    def csv_headers(self, separator: str = ",") -> str:
        return separator.join(
            [
                "domain_alignment_id",
                "feature_id",
            ]
        )

    @classmethod
    def from_query_result(
        cls,
        query_result: HmmscanQueryResult,
    ) -> List["IntersectionMapTableRecord"]:
        if not query_result.intersections:
            return []

        records = []

        for intersection in query_result.intersections:
            records.append(
                cls(
                    domain_alignment_id=intersection.domain_alignment_id,
                    feature_id=intersection.feature_id,
                )
            )

        return records


@dataclass
class GroupTableRecord:
    group_id: str
    group_start_position: int
    group_end_position: int
    domain_alignments_in_group: List[str]
    domain_alignments_in_pos_strand: List[str]
    domain_alignments_in_neg_strand: List[str]

    def __post_init__(self):
        self.n_domain_alignments = len(self.domain_alignments_in_group)
        self.n_domain_alignments_in_pos_strand = len(
            self.domain_alignments_in_pos_strand
        )
        self.n_domain_alignments_in_neg_strand = len(
            self.domain_alignments_in_neg_strand
        )

    def to_csv(self, separator: str = ",") -> str:
        return separator.join(
            [
                self.group_id,
                str(self.group_start_position),
                str(self.group_end_position),
                str(self.n_domain_alignments),
                str(self.n_domain_alignments_in_pos_strand),
                str(self.n_domain_alignments_in_neg_strand),
                ";".join(self.domain_alignments_in_group),
                ";".join(self.domain_alignments_in_pos_strand),
                ";".join(self.domain_alignments_in_neg_strand),
            ]
        )

    def csv_headers(self, separator: str = ",") -> str:
        return separator.join(
            [
                "group_id",
                "group_start_position",
                "group_end_position",
                "n_domain_alignments",
                "n_domain_alignments_in_pos_strand",
                "n_domain_alignments_in_neg_strand",
                "domain_alignments_in_group",
                "domain_alignments_in_pos_strand",
                "domain_alignments_in_neg_strand",
            ]
        )

    @classmethod
    def from_query_result(
        cls,
        query_result: HmmscanQueryResult,
    ) -> "GroupTableRecord | None":
        if not query_result.group:
            return None

        return cls(
            group_id=query_result.group.id,
            group_start_position=query_result.group.start,
            group_end_position=query_result.group.end,
            domain_alignments_in_group=query_result.group.domain_alignments_in_group,
            domain_alignments_in_pos_strand=query_result.group.domain_alignments_in_pos_strand,
            domain_alignments_in_neg_strand=query_result.group.domain_alignments_in_neg_strand,
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

    # Optional arguments
    parser.add_argument(
        "--kind",
        type=str,
        default="sequence",
        choices=["sequence", "domain", "gff", "map", "group"],
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

    # Parse the serialized query results from JSON to dataclasses again.
    query_results = []
    for line in file_handle:
        try:
            query_result = HmmscanQueryResult.from_json(json.loads(line))
        except UnexpectedQueryIdFormat as e:
            logger.error(e)
            sys.exit(1)

        query_results.append(query_result)

    # Make a domain set to be used in the binary representation of the domain hits
    # in the sequence-centric output.
    domain_set = set()
    for query_result in query_results:
        for domain_hit in query_result.domain_hits:
            domain_set.add(domain_hit.name)

    table_records = []
    match config.kind:
        case "sequence":
            table_records = [
                SequenceCentricTableRecord.from_query_result(query_result, domain_set)
                for query_result in query_results
            ]

        case "domain":
            for query_result in query_results:
                table_records.extend(
                    DomainAlignmentCentricTableRecord.from_domain_hit(query_result)
                )

        case "gff":
            for query_result in query_results:
                table_records.extend(
                    IntersectingGffFeatureTableRecord.from_query_result(query_result)
                )

        case "map":
            for query_result in query_results:
                table_records.extend(
                    IntersectionMapTableRecord.from_query_result(query_result)
                )

        case "group":
            if not any(query_result.group for query_result in query_results):
                logger.error("No groups found in the query results.")
                sys.exit(1)

            table_records = [
                GroupTableRecord.from_query_result(query_result)
                for query_result in query_results
                if query_result.group
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
