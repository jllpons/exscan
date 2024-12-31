#!/usr/bin/env python

"""
domtblop parser: Parse and serialize hmmscan domtblout query results to JSON.

Usage: domtblop.py parse <domtblout> --seq-type <sequence type>

Arguments:
    domtblout   Path to the hmmscan domtblout file.
    --seq-type  Type of the sequence used in the hmmscan search.
                Sequence types: {dna, rna, protein}


Options:
    -l --loglevel   Log level (default: INFO).
    -h --help       Show this screen.
    -v --version    Show version.
"""

import argparse
from dataclasses import dataclass, field
from datetime import datetime
import json
import logging
import os
import sys
from typing import (
    Dict,
    Generator,
    List,
    Tuple,
)

from Bio.SearchIO import parse
from Bio.SearchIO._model import QueryResult

from domtblop import __version__
from domtblop_utils import (
    setup_logger,
)


SCHEMA_VERSION = "0.1.0"


class UnexpectedQueryIdFormat(Exception):
    """
    Raised when the hit id is not in the expected format.
    """

    pass


class CustomEncoder(json.JSONEncoder):
    """
    A workaround to serialize dataclasses to JSON if they have a `to_json` method.
    """

    def default(self, obj):
        if hasattr(obj, "to_json"):
            return obj.to_json()
        return super().default(obj)


@dataclass
class QueryResultMetadata:
    """
    Holds metadata information about a query result.

    Attributes:
        schema_version (str): Version of the schema used to serialize the query result.
    """

    schema_version: str = SCHEMA_VERSION
    domtblop_version: str = __version__
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())
    last_modified_at: str = field(default_factory=lambda: datetime.now().isoformat())
    parameters_used: list = field(default_factory=list)

    def to_json(self) -> Dict:
        """
        Serialize the metadata into a JSON struct.

        Returns:
            Dict: JSON struct of the metadata.
        """
        return {
            "schema_version": self.schema_version,
            "domtblop_version": self.domtblop_version,
            "created_at": self.created_at,
            "last_modified_at": self.last_modified_at,
            "parameters_used": [i for i in self.parameters_used],
        }

    @classmethod
    def from_json(cls, json_data: Dict) -> "QueryResultMetadata":
        """
        Deserialize a JSON struct into a QueryResultMetadata.

        Args:
            json_data (Dict): JSON string to deserialize.

        Returns:
            QueryResultMetadata: Deserialized QueryResultMetadata.

        Raises:
            KeyError: If a required key is missing in the JSON struct.
        """

        try:
            schema_version = json_data["schema_version"]
            domtblop_version = json_data["domtblop_version"]
            created_at = json_data["created_at"]
            last_modified_at = json_data["last_modified_at"]
            parameters_used = json_data["parameters_used"]

        except KeyError as e:
            raise KeyError(f"Missing key: {e}")

        metadata = cls(
            schema_version=schema_version,
            domtblop_version=domtblop_version,
            created_at=created_at,
            last_modified_at=last_modified_at,
            parameters_used=parameters_used,
        )

        return metadata


@dataclass
class ParentSequence:
    """
    Holds information about the parent sequence of a hit.

    Attributes:
        sequence_id (str): Identifier of the parent sequence.
        start (int): Start position of the hit in the parent sequence.
        end (int): End position of the hit in the parent sequence.
        strand (str): Strand information ('+' or '-').
        frame (int): Reading frame of the hit.
    """

    sequence_id: str
    start: int
    end: int
    strand: str
    frame: int

    def to_json(self) -> Dict:
        """
        Serialize the parent sequence into a JSON struct.

        Returns:
            Dict: JSON struct of the parent sequence.
        """
        return {
            "sequence_id": self.sequence_id,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "frame": self.frame,
        }

    @classmethod
    def from_json(cls, json_data: Dict) -> "ParentSequence":
        """
        Deserialize a JSON struct into a ParentSequence.

        Args:
            json_data (Dict): JSON string to deserialize.

        Returns:
            ParentSequence: Deserialized ParentSequence.

        Raises:
            KeyError: If a required key is missing in the JSON struct.
        """

        try:
            sequence_id = json_data["sequence_id"]
            start = json_data["start"]
            end = json_data["end"]
            strand = json_data["strand"]
            frame = json_data["frame"]

        except KeyError as e:
            raise KeyError(f"Missing key: {e}")

        parent_sequence = cls(
            sequence_id=sequence_id,
            start=start,
            end=end,
            strand=strand,
            frame=frame,
        )

        return parent_sequence


@dataclass
class GffFeature:
    """
    Represents a GFF feature that intersects with a domain alignment fragment.

    Attributes:
        seqid (str): Sequence identifier.
        source (str): Source of the feature.
        type_ (str): Type of the feature.
        start (int): Start position of the feature.
        end (int): End position of the feature.
        score (float | str): Score of the feature.
        strand (str): Strand of the feature.
        phase (int | str): Phase of the feature.
        attributes (Dict[str, str]): Additional attributes of the feature
    """

    seqid: str
    source: str
    type_: str
    start: int
    end: int
    score: float | str
    strand: str
    phase: int | str
    attributes: Dict[str, str]

    def to_json(self) -> Dict:
        """
        Serialize the GFF feature into a JSON struct.

        Returns:
            Dict: JSON struct of the GFF feature.
        """
        return {
            "seqid": self.seqid,
            "source": self.source,
            "type": self.type_,
            "start": self.start,
            "end": self.end,
            "score": self.score,
            "strand": self.strand,
            "phase": self.phase,
            "attributes": self.attributes,
        }

    @classmethod
    def from_json(cls, json_data: Dict) -> "GffFeature":
        """
        Deserialize a JSON struct into a GffIntersectingFeature.

        Args:
            json_data (Dict): JSON string to deserialize.

        Returns:
            GffIntersectingFeature: Deserialized GffIntersectingFeature.

        Raises:
            KeyError: If a required key is missing in the JSON struct.
        """

        try:
            seqid = json_data["seqid"]
            source = json_data["source"]
            type_ = json_data["type"]
            start = json_data["start"]
            end = json_data["end"]
            score = json_data["score"]
            strand = json_data["strand"]
            phase = json_data["phase"]
            attributes = json_data["attributes"]

            if seqid.startswith('"') and seqid.endswith('"'):
                seqid = seqid[1:-1]

        except KeyError as e:
            raise KeyError(f"Missing key: {e}")

        feature = cls(
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

        return feature

    @classmethod
    def from_gff3(cls, line: str) -> "GffFeature":
        """
        Create a Gff3Feature instance from a GFF3 line.

        Args:
            line (str): A GFF3 line.

        Returns:
            Gff3Feature: A Gff3Feature instance.
        """

        def _parse_attributes(attributes: str) -> dict:
            fields = attributes.split(";")
            return {
                field.split("=")[0].replace('"', ""): "=".join(field.split("=")[1:])
                for field in fields
                if field
            }

        fields = line.strip().split("\t")
        attributes = _parse_attributes(fields[8])

        return cls(
            seqid=fields[0],
            source=fields[1],
            type_=fields[2],
            start=int(fields[3]),
            end=int(fields[4]),
            score=float(fields[5]) if fields[5] != "." else ".",
            strand=fields[6],
            phase=int(fields[7]) if fields[7] != "." else ".",
            attributes=attributes,
        )

    def to_gff3(self) -> str:
        """
        Serialize the GFF feature into a GFF3 line.

        Returns:
            str: GFF3 line of the GFF feature.
        """
        attributes = ";".join(f"{k}={v}" for k, v in self.attributes.items())

        return "\t".join(
            [
                self.seqid,
                self.source,
                self.type_,
                str(self.start),
                str(self.end),
                str(self.score),
                self.strand,
                str(self.phase),
                attributes,
            ]
        )


def print_gff3_feature(feature: GffFeature) -> None:
    """
    Print a Gff3Feature object in GFF3 format.

    Args:
        feature (Gff3Feature): A Gff3Feature object.
    """
    print(
        "\t".join(
            [
                feature.seqid,
                feature.source,
                feature.type_,
                str(feature.start),
                str(feature.end),
                str(feature.score),
                feature.strand,
                str(feature.phase),
                ";".join(f"{k}={v}" for k, v in feature.attributes.items()),
            ]
        )
    )


@dataclass
class DomainAlignmentFragment:
    """
    Represents a fragment of a domain alignment from hmmscan results.

    Attributes:
        domain_start (int): Start position in the domain model.
        domain_end (int): End position in the domain model.
        domain_strand (str): Strand in the domain model (if applicable).
        sequence_start (int): Start position in the query sequence.
        sequence_end (int): End position in the query sequence.
        sequence_strand (str): Strand in the query sequence (if applicable).
        sequence_start_in_parent (None | int): Start position in the parent sequence.
        sequence_end_in_parent (None | int): End position in the parent sequence
        sequence (None | str): The sequence of the fragment.
    """

    domain_start: int
    domain_end: int
    domain_strand: str
    sequence_start: int
    sequence_end: int
    sequence_strand: str

    sequence_start_in_parent: None | int = None
    sequence_end_in_parent: None | int = None
    sequence: None | str = None

    def to_json(self) -> Dict:
        """
        Serialize the hmmscan HSP fragment into a JSON struct.

        Returns:
            Dict: JSON struct of the hmmscan HSP fragment.
        """
        return {
            "domain_start": self.domain_start,
            "domain_end": self.domain_end,
            "domain_strand": self.domain_strand,
            "sequence_start": self.sequence_start,
            "sequence_end": self.sequence_end,
            "sequence_strand": self.sequence_strand,
            "sequence_start_in_parent": self.sequence_start_in_parent,
            "sequence_end_in_parent": self.sequence_end_in_parent,
            "sequence": self.sequence,
        }

    @classmethod
    def from_json(cls, json_data: Dict) -> "DomainAlignmentFragment":
        """
        Deserialize a JSON struct into a HmmscanHspFragment.

        Args:
            json_data (Dict): JSON string to deserialize.

        Returns:
            DomainAlignmentFragment: Deserialized DomainAlignmentFragment.

        Raises:
            KeyError: If a required key is missing in the JSON struct.
        """

        try:
            domain_start = json_data["domain_start"]
            domain_end = json_data["domain_end"]
            domain_strand = json_data["domain_strand"]
            sequence_start = json_data["sequence_start"]
            sequence_end = json_data["sequence_end"]
            sequence_strand = json_data["sequence_strand"]
            sequence_start_in_parent = json_data["sequence_start_in_parent"]
            sequence_end_in_parent = json_data["sequence_end_in_parent"]
            sequence = json_data["sequence"]

        except KeyError as e:
            raise KeyError(f"Missing key: {e}")

        hspfrag = cls(
            domain_start=domain_start,
            domain_end=domain_end,
            domain_strand=domain_strand,
            sequence_start=sequence_start,
            sequence_end=sequence_end,
            sequence_strand=sequence_strand,
            sequence_start_in_parent=sequence_start_in_parent,
            sequence_end_in_parent=sequence_end_in_parent,
            sequence=sequence,
        )

        return hspfrag


@dataclass
class HmmscanDomainAlignment:
    """
    Represents an alignment of a domain within a query sequence.

    Attributes:
        independent_evalue (float): i-Evalue for the domain alignment.
        conditional_evalue (float): c-Evalue for the domain alignment.
        average_alignment_accuracy (float): Average alignment accuracy.
        bias (float): Bias score for the alignment.
        bit_score (float): Bit score for the domain alignment.
        domain_number (int): The index of the domain in the query sequence.
        envelope_start (int): Start of the envelope region in the sequence.
        envelope_end (int): End of the envelope region in the sequence.
        alignment_fragments (List[DomainAlignmentFragment]): List of alignment fragments.
    """

    independent_evalue: float
    conditional_evalue: float
    average_alignment_accuracy: float
    bias: float
    bit_score: float
    domain_number: int
    envelope_start: int
    envelope_end: int
    alignment_fragments: List[DomainAlignmentFragment] = field(default_factory=list)

    def to_json(self) -> Dict:
        """
        Serialize the hmmscan HSP into a JSON struct.

        Returns:
            Dict: JSON struct of the hmmscan HSP.
        """
        return {
            "independent_evalue": self.independent_evalue,
            "conditional_evalue": self.conditional_evalue,
            "average_alignment_accuracy": self.average_alignment_accuracy,
            "bias": self.bias,
            "bit_score": self.bit_score,
            "domain_number": self.domain_number,
            "envelope_start": self.envelope_start,
            "envelope_end": self.envelope_end,
            "alignment_fragments": [i.to_json() for i in self.alignment_fragments],
        }

    @classmethod
    def from_json(cls, json_data: Dict) -> "HmmscanDomainAlignment":
        """
        Deserialize a JSON struct into a HmmscanHsp.

        Args:
            json_data (Dict): JSON string to deserialize.

        Returns:
            HmmscanDomainAlignment: Deserialized HmmscanDomainAlignment.

        Raises:
            KeyError: If a required key is missing in the JSON struct.
        """

        try:
            independent_evalue = json_data["independent_evalue"]
            conditional_evalue = json_data["conditional_evalue"]
            average_alignment_accuracy = json_data["average_alignment_accuracy"]
            bias = json_data["bias"]
            bit_score = json_data["bit_score"]
            domain_number = json_data["domain_number"]
            envelope_start = json_data["envelope_start"]
            envelope_end = json_data["envelope_end"]
            alignment_fragments = json_data["alignment_fragments"]

        except KeyError as e:
            raise KeyError(f"Missing key: {e}")

        domain_alignment = cls(
            independent_evalue=independent_evalue,
            conditional_evalue=conditional_evalue,
            average_alignment_accuracy=average_alignment_accuracy,
            bias=bias,
            bit_score=bit_score,
            domain_number=domain_number,
            envelope_start=envelope_start,
            envelope_end=envelope_end,
            alignment_fragments=[
                DomainAlignmentFragment.from_json(i) for i in alignment_fragments
            ],
        )

        return domain_alignment

    def has_domain_alignment_fragments(self) -> bool:
        """
        Check if the domain alignment has fragments.

        Returns:
            bool: True if the domain alignment has fragments, False otherwise.
        """
        return any(self.alignment_fragments)


@dataclass
class HmmscanDomainHit:
    """
    Represents a domain hit from an hmmscan result.

        Attributes:
            accession (str): Accession number of the domain.
            name (str): Name or identifier of the domain.
            description (str): Description of the domain.
            full_sequence_evalue (float): E-value for the full sequence.
            full_sequence_score (float): Score for the full sequence.
            bias (float): Bias score for the domain hit.
            domain_alignments (List[HmmscanDomainAlignment]): List of domain alignments.
    """

    accession: str
    name: str
    description: str
    full_sequence_evalue: float
    full_sequence_score: float
    bias: float
    domain_alignments: List[HmmscanDomainAlignment] = field(default_factory=list)

    def to_json(self) -> Dict:
        """
        Serialize the hmmscan into a JSON struct.

        Returns:
            Dict: JSON struct of the hmmscan.
        """
        return {
            "accession": self.accession,
            "name": self.name,
            "description": self.description,
            "full_sequence_evalue": self.full_sequence_evalue,
            "full_sequence_score": self.full_sequence_score,
            "bias": self.bias,
            "domain_alignments": [i.to_json() for i in self.domain_alignments],
        }

    @classmethod
    def from_json(cls, json_data: Dict) -> "HmmscanDomainHit":
        """
        Deserialize a JSON struct into a HmmscanDomainHit.

        Args:
            json_data (Dict): JSON string to deserialize.

        Returns:
            HmmscanDomainHit: Deserialized HmmscanDomainHit.

        Raises:
            KeyError: If a required key is missing in the JSON struct.
        """

        try:
            accession = json_data["accession"]
            name = json_data["name"]
            description = json_data["description"]
            full_sequence_evalue = json_data["full_sequence_evalue"]
            full_sequence_score = json_data["full_sequence_score"]
            bias = json_data["bias"]
            domain_alignments = json_data["domain_alignments"]

        except KeyError as e:
            raise KeyError(f"Missing key: {e}")

        hit = cls(
            accession=accession,
            name=name,
            description=description,
            full_sequence_evalue=full_sequence_evalue,
            full_sequence_score=full_sequence_score,
            bias=bias,
            domain_alignments=[
                HmmscanDomainAlignment.from_json(i) for i in domain_alignments
            ],
        )

        return hit

    def has_domain_alignments(self) -> bool:
        """
        Check if the domain hit has alignments.

        Returns:
            bool: True if the domain hit has alignments, False otherwise.
        """
        return any(self.domain_alignments)


@dataclass
class Group:
    """
    Holds information about a group of HmmscanQueryResults that are found within
    *n* distance form each other according to the parent sequence start and end
    coordinates

    Atributes:
        group_id (str):
        group_start (int):
        group_end (int):
        n_hits_pos_strand (int):
        n_hits_neg_strand (int):
    """

    id: str
    start: int
    end: int
    n_hits_on_pos_strand: int
    n_hits_on_neg_strand: int

    def to_json(self) -> Dict:
        """
        Serialize the hit into a JSON struct.

        Returns:
            Dict: JSON struct of the hit.
        """

        return {
            "id": self.id,
            "start": self.start,
            "end": self.end,
            "n_hits_on_pos_strand": self.n_hits_on_pos_strand,
            "n_hits_on_neg_strand": self.n_hits_on_neg_strand,
        }

    @classmethod
    def from_json(cls, json_data: Dict) -> "Group":
        """
        Deserialize a JSON struct into a Group.

        Args:
            json_data (Dict): JSON string to deserialize.

        Returns:
            GroupedQueryResults: Deserialized GroupedQueryResults.

        Raises:
            KeyError: If a required key is missing in the JSON struct.
        """

        try:
            id = json_data["id"]
            group_start = json_data["start"]
            group_end = json_data["end"]
            n_hits_on_pos_strand = json_data["n_hits_on_pos_strand"]
            n_hits_on_neg_strand = json_data["n_hits_on_neg_strand"]

        except KeyError as e:
            raise KeyError(f"Missing key: {e}")

        grouped_query_result = cls(
            id=id,
            start=group_start,
            end=group_end,
            n_hits_on_pos_strand=n_hits_on_pos_strand,
            n_hits_on_neg_strand=n_hits_on_neg_strand,
        )

        return grouped_query_result


@dataclass
class HmmscanQueryResult:
    """
    Represents the result of an hmmscan query for a single sequence.

    Attributes:
        query_id (str): Identifier of the query sequence.
        aminoacid_sequence (str): Amino acid sequence of the query.
        nucleotide_sequence (str): Nucleotide sequence of the query.
        parent_sequence (ParentSequence): Information about the parent sequence.
        domain_hits (List[HmmscanDomainHit]): List of domain hits found in the query.
        group (None | GroupedQueryResults): Data about a set of hits within certain distance.
    """

    metadata: QueryResultMetadata
    query_id: str
    sequence_type: str
    protein_sequence: str | None = None
    source_sequence: str | None = None
    parent_sequence: ParentSequence | None = None
    domain_hits: List[HmmscanDomainHit] = field(default_factory=list)

    intersecting_features: None | List[GffFeature] = None
    group: None | Group = None

    def __post_init__(self):
        """
        The code within this method is executed after the dataclass is
        initialized. Parsing the parent sequence information from the hit id.

        If the hit has already been provided with a parent sequence, the
        method will not attempt to parse the parent sequence from the hit id.

        If the hit id is not in the expected format, an UnexpectedHitIdFormat

        Returns:
            None

        Raises:
            UnexpectedHitIdFormat: If the hit id is not in the expected format

        Notes:
            In order to know the location of a hit within a parent sequence,
            the fasta header of the parent sequence is expected to have the
            following format:
                <parent_id>;locOrf=<start>..<end>
            where:
                parent_id: The id of the parent sequence.
                start: The start position of the hit.
                end: The end position of the hit.
        """
        self.check_sequence_type()

        if self.sequence_type in ["dna", "rna"]:
            if "_frame=" not in self.query_id:
                raise UnexpectedQueryIdFormat(
                    "Query id must contain '_frame=', '_start=', and '_end='"
                )

            parts = self.query_id.split("_")
            sequence_id = "_".join(parts[:-3])

            frame_part, start_part, end_part = parts[-3:]
            try:
                frame = int(frame_part.split("=")[1])
                start = int(start_part.split("=")[1])
                end = int(end_part.split("=")[1])
            except ValueError:
                raise UnexpectedQueryIdFormat(
                    "Start and end positions must be integers"
                )

            strand = "+" if frame >= 0 else "-"

            self.parent_sequence = ParentSequence(
                sequence_id=sequence_id,
                start=start,
                end=end,
                strand=strand,
                frame=frame,
            )

    def to_json(self) -> Dict:
        """
        Serialize the hit into a JSON struct.

        Returns:
            Dict: JSON struct of the hit.
        """

        return {
            "metadata": self.metadata.to_json(),
            "query_id": self.query_id,
            "sequence_type": self.sequence_type,
            "protein_sequence": self.protein_sequence,
            "source_sequence": self.source_sequence,
            "parent_sequence": self.parent_sequence.to_json()
            if self.parent_sequence
            else None,
            "domain_hits": [i.to_json() for i in self.domain_hits],
            "intersecting_features": [i.to_json() for i in self.intersecting_features]
            if self.intersecting_features
            else None,
            "group": self.group.to_json() if self.group else None,
        }

    @classmethod
    def from_json(cls, json_data: Dict) -> "HmmscanQueryResult":
        """
        Deserialize a JSON struct into a HmmscanQueryResult.

        Args:
            json_data (Dict): JSON string to deserialize.

        Returns:
            HmmscanQueryResult: Deserialized HmmscanQueryResult.

        Raises:
            KeyError: If a required key is missing in the JSON struct.

        Example:

            # Serialize a Hit instance to JSON
            query_result = HmmscanQueryResult(
                query_id="ABC123;locOrf=100..200",
                #...
            )
            json_data = query_result.to_json()

            print(json.dumps(json_data))
            # >> {"query_id": "ABC123;locOrf=100..200", ...}

            # Deserialize from JSON to a HmmscanQueryResult instance
            query_result = HmmscanQueryResult.from_json(json_data)
        """

        try:
            metadata = json_data["metadata"]
            query_id = json_data["query_id"]
            sequence_type = json_data["sequence_type"]
            protein_sequence = json_data["protein_sequence"]
            source_sequence = json_data["source_sequence"]
            domain_hits = json_data["domain_hits"]
            intersecting_features = json_data["intersecting_features"]
            group = json_data["group"]

        except KeyError as e:
            raise KeyError(f"Missing key: {e}")

        query_result = cls(
            metadata=QueryResultMetadata.from_json(metadata),
            query_id=query_id,
            sequence_type=sequence_type,
            protein_sequence=protein_sequence,
            source_sequence=source_sequence,
            domain_hits=[HmmscanDomainHit.from_json(i) for i in domain_hits],
            intersecting_features=[
                GffFeature.from_json(i) for i in intersecting_features
            ]
            if intersecting_features
            else None,
            group=Group.from_json(group) if group else None,
        )

        return query_result

    def check_sequence_type(self) -> None:
        """
        Check if the sequence type is valid.

        Returns:
            None
        """
        if self.sequence_type not in ["dna", "rna", "protein"]:
            raise ValueError("Invalid sequence type: {}".format(self.sequence_type))

    def update_modified_at(self) -> None:
        """
        Update the modified_at field of the metadata.

        Returns:
            None
        """
        self.metadata.last_modified_at = datetime.now().isoformat()

    def has_domain_hits(self) -> bool:
        """
        Check if the query result has domain hits.

        Returns:
            bool: True if the query result has domain hits, False otherwise.
        """
        return any(self.domain_hits)

    def has_intersecting_gff_features(self) -> bool:
        """
        Check if the query result has intersecting GFF features.

        Returns:
            bool: True if the query result has intersecting GFF features, False otherwise.
        """
        if not self.intersecting_features:
            return False

        return any(self.intersecting_features)

    def compute_domain_aligment_fragment_absolut_positions(self) -> None:
        """
        Compute the absolute positions of the domain alignments in the parent sequence.

        Returns:
            None
        """
        if self.sequence_type not in ["dna", "rna"]:
            return

        strand = self.parent_sequence.strand

        for domain_hit in self.domain_hits:
            for domain_alignment in domain_hit.domain_alignments:
                for fragment in domain_alignment.alignment_fragments:
                    # Fragment positions are amino acid based
                    if strand == "+":
                        fragment.sequence_start_in_parent = (
                            self.parent_sequence.start + (fragment.sequence_start * 3)
                        )
                        fragment.sequence_end_in_parent = self.parent_sequence.start + (
                            fragment.sequence_end * 3
                        )

                    elif strand == "-":
                        # feature_lenght = (fragment.sequence_end - fragment.sequence_start) * 3
                        fragment.sequence_start_in_parent = self.parent_sequence.end - (
                            fragment.sequence_end * 3
                        )
                        fragment.sequence_end_in_parent = self.parent_sequence.end - (
                            fragment.sequence_start * 3
                        )

    def mk_domain_alignment_fragment_sequences(self) -> None:
        """
        Make the sequences for each domain alignment fragment.

        Args:
            sequence (str): The parent sequence.

        Returns:
            None
        """
        if not self.protein_sequence:
            raise ValueError("Missing protein sequence")
        for domain_hit in self.domain_hits:
            for domain_alignment in domain_hit.domain_alignments:
                for fragment in domain_alignment.alignment_fragments:
                    fragment.sequence = self.protein_sequence[
                        fragment.sequence_start : fragment.sequence_end
                    ]


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
        "domtblout",
        metavar="<domtblout>",
        type=str,
        nargs="?",
    )
    parser.add_argument(
        "--seq-type",
        metavar="<sequence type>",
        type=str,
        choices=["dna", "rna", "protein"],
        required=True,
    )

    # Optional arguments
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

    if not config.domtblout:
        parser.print_help()
        logger.error("Missig required argument <domtblout>.")
        raise ValueError

    if not os.path.exists(config.domtblout):
        logger.error("File not found: {}".format(config.domtblout))
        raise FileNotFoundError

    if not config.seq_type:
        parser.print_help()
        logger.error("Missig required argument <seq-type>.")
        raise

    if config.seq_type not in ["dna", "rna", "protein"]:
        parser.print_help()
        logger.error("Invalid sequence type: {}".format(config.seq_type))
        raise ValueError

    return config, logger


def parse_hmmscan_domtblout(hmmscan_domtblout: str) -> Generator:
    """
    Parse hmmscan domtblout file. Return a generator of query results.

    Args:
        hmmscan_domtblout (str): Path to hmmscan domtblout file.

    Returns:
                   Generator of query results.
                   Each query result is a Bio.SearchIO object
                   For more information, see:
                     <https://biopython.org/docs/1.75/api/Bio.SearchIO.html>

    Raises:
        AssertionError: If the file does not exist.
    """

    if not os.path.isfile(hmmscan_domtblout):
        raise AssertionError("File does not exist: {}".format(hmmscan_domtblout))

    queryresult_generator = parse(hmmscan_domtblout, "hmmscan3-domtab")

    return queryresult_generator


def parse_query_result(
    query_result: QueryResult,
    seq_type: str,
) -> HmmscanQueryResult:
    """
    Parse a QueryResult object from a hmmscan domtbl file.

    Args:
        query_result (QueryResult): A QueryResult object from a hmmscan domtbl
        seq_type (str): The type of the sequence used in the hmmscan search.

    Returns:
        HmmscanQueryResult: A parsed HmmscanQueryResult object.

    Raises:
        UnexpectedHitIdFormat: If the hit id is not in the expected format.
    """
    if seq_type not in ["dna", "rna", "protein"]:
        raise ValueError("Invalid sequence type: {}".format(seq_type))

    query_id = query_result.id

    domain_hits = []
    for domain_hit in query_result:
        accession = domain_hit.accession
        name = domain_hit.id
        description = domain_hit.description
        full_sequence_evalue = domain_hit.evalue
        full_sequence_score = domain_hit.bitscore
        bias = domain_hit.bias
        domain_alignments = []

        for alignment in domain_hit:
            independent_evalue = alignment.evalue
            conditional_evalue = alignment.evalue_cond
            average_aligment_accuracy = alignment.acc_avg
            bias = alignment.bias
            bit_score = alignment.bitscore
            domain_number = alignment.domain_index
            envelop_start = alignment.env_start
            envelope_end = alignment.env_end
            alignment_fragments = []

            for fragment in alignment:
                domain_start = fragment.hit_start
                domain_end = fragment.hit_end
                domain_strand = fragment.hit_strand
                sequence_start = fragment.query_start
                sequence_end = fragment.query_end
                sequence_strand = fragment.query_strand

                alignment_fragments.append(
                    DomainAlignmentFragment(
                        domain_start=domain_start,
                        domain_end=domain_end,
                        domain_strand=domain_strand,
                        sequence_start=sequence_start,
                        sequence_end=sequence_end,
                        sequence_strand=sequence_strand,
                    )
                )

            domain_alignments.append(
                HmmscanDomainAlignment(
                    independent_evalue=independent_evalue,
                    conditional_evalue=conditional_evalue,
                    average_alignment_accuracy=average_aligment_accuracy,
                    bias=bias,
                    bit_score=bit_score,
                    domain_number=domain_number,
                    envelope_start=envelop_start,
                    envelope_end=envelope_end,
                    alignment_fragments=alignment_fragments,
                )
            )

        domain_hits.append(
            HmmscanDomainHit(
                accession=accession,
                name=name,
                description=description,
                full_sequence_evalue=full_sequence_evalue,
                full_sequence_score=full_sequence_score,
                bias=bias,
                domain_alignments=domain_alignments,
            )
        )

    metadata = QueryResultMetadata()

    parsed_query_result = HmmscanQueryResult(
        metadata=metadata,
        query_id=query_id,
        sequence_type=seq_type,
        domain_hits=domain_hits,
    )

    # Call any additional methods to compute or modify the data
    parsed_query_result.compute_domain_aligment_fragment_absolut_positions()

    return parsed_query_result


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
        queryresult_generator = parse_hmmscan_domtblout(config.domtblout)
    except AssertionError as e:
        logger.error(f"{e}")
        sys.exit(1)

    for query_result in queryresult_generator:
        logger.debug(f"Parsing query result: {query_result.id}")

        try:
            parsed_query_result = parse_query_result(query_result, config.seq_type)
        except (UnexpectedQueryIdFormat, ValueError) as e:
            logger.error("Error: {}".format(e))
            exit(1)

        logger.debug(f"Succesfully parsed query result: {query_result.id}")

        try:
            print(json.dumps(parsed_query_result, cls=CustomEncoder))

        # <https://docs.python.org/3/library/signal.html#note-on-sigpipe>
        except BrokenPipeError:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())
            sys.exit(1)
