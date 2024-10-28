#!/usr/local/bin/gawk -f

# `getorg` from EMBOSS package outpus like this:
# ```
# >myFavoriteGene [1 - 100]
# KAPPSKRKTRNNHDVRLGLIWAGPIPEELEDLTELKELWLDNNNLTGKGCWIRHDVSSTS
# ```
# ```
#
# It can also contain more information like:
# ```
# >LEXA_ECOLI_1 [11 - 52] LexA repressor OS=Escherichia coli (strain K12) OX=83333 GN=lexA PE=1 SV=1 to a 606 base sequence of most likely codons.
# KAPPSKRKTRNNHDVRLGLIWAGPIPEELEDLTELKELWLDNNNLTGKGCWIRHDVSSTS
# ```
#
#  This script should convert the first to:
#  ```
#  >myFavoriteGene;locOrf=1..100
#  ```
#  And the second to:
#  ```
#  >LEXA_ECOLI_1;locOrf=11..52 LexA repressor OS=Escherichia coli (strain K12) OX=83333 GN=lexA PE=1 SV=1
#  ```
#
#  Main motivation is to conserve the cordinates in the hmmscan --domtblout output
{
    if(!/^>/){
        print $0;
    }
    else if(/^>/){
        # The id is the first word after the '>'
        id = $1;
        # `[N`
        int1 = $2;
        hyphen = $3;
        # `M]`
        int2 = $4;
        # Remove the brackets
        gsub(/\[|\]/, "", int1);
        gsub(/\[|\]/, "", int2);

        printf "%s;locOrf=%s..%s\n", id, int1, int2;
    }
}
