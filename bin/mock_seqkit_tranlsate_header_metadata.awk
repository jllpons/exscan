#!/usr/bin/awk -f

# When a nulceotide sequence is translated using the `seqkit` cli tool
# with the `transalte` command, the header of the sequence is modified
# to include the translation frame and the coordinates of the nucleotide sequence.
#
# I.e.:
#
# From:
# ```fasta
# >myFavoriteTranslatedGene Some other metadata;GO:0000000;PFAM:PF00000...
# KAPPSKRKTRNNHDVRLGLIWAGPIPEELEDLTELKELWLDNNNLTGKGCWIRHDVSSTS
# ```
#
# To:
# ```fasta
# >myFavoriteTranslatedGene_frame=1_start=33_end=99
# KAPPSKRKTRNNHDVRLGLIWAGPIPEELEDLTELKELWLDNNNLTGKGCWIRHDVSSTS
# ```
#
# Some downstream processes of this pipeline rely on this information.
# When the pipeline is used for already translated sequences, or mRNAs,
# the header metadata is "replicated" using this script.


{
    if ($0 ~ /^>/) {
        # If we have a previous sequence, process it
        if (sequence != "") {
            # Print the header with metadata
            print header "_frame=" frame "_start=" start "_end=" length(sequence)
            # Print the sequence
            print sequence
        }
        # Update header and reset sequence
        header = $0
        split(header, header_parts, " ")
        header = header_parts[1]
        frame = "0"   # Indicate that the sequence is already translated
        start = "0"
        sequence = ""
    } else {
        # Append the current line to the sequence
        sequence = sequence $0
    }
}

END {
    # Print the last sequence if it exists
    if (sequence != "") {
        print header "_frame=" frame "_start=" start "_end=" length(sequence)
        print sequence
    }
}
