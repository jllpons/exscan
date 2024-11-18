#!/usr/bin/env python

"""
domtblop: Perform operations on query results from a hmmscan domtblout files.

Usage: domtblop.py <operation> [options]

    Parsing and Serialization:
        parser   Parse and serialize hmmscan domtblout query results into JSON structs

    Operations on Serialized Domtblout Query Results:
        filter     Apply filters to serialized structs (e.g. E-value, score)
        group      Group domain hits found within a certain distance of each other.
        gffcmp     Compare domain hits or groups to a GFF file.
        bamcmp     Compare domain hits or groups to a BAM file.
        togff      Convert serialized domain hits or groups to GFF format.
        tocsv      Convert serialized domain hits or groups to CSV format.

    Utilities:
        addseq     Add the amino acid sequence to each query result.
        dummy      Print some serialized dummy query results for testing purposes.
        help       Show this help message and exit.

Options:
    -h, --help   Show this help message and exit.
"""


import sys


def main():

    argn = len(sys.argv)
    if argn == 1:
        print(__doc__)
        sys.exit(1)

    cmd = sys.argv[1]
    cwd = __file__.rsplit("/", 1)[0]

    if cmd == "-h" or cmd == "--help" or cmd == "help":
        print(__doc__)
        sys.exit(1)

    elif cmd == "parser":
        from domtblop_parser import run
        run(sys.argv[2:])

    elif cmd == "filter":
        from domtblop_filter import run
        run(sys.argv[2:])

    elif cmd == "group":
        from domtblop_group import run
        run(sys.argv[2:])

    elif cmd == "gffcmp":
        subcmd = [f"{cwd}/domtblop_gffcmp.py"] + sys.argv[2:]

    elif cmd == "bamcmp":
        subcmd = [f"{cwd}/domtblop_bamcmp.py"] + sys.argv[2:]

    elif cmd == "togff":
        from domtblop_togff import run
        run(sys.argv[2:])

    elif cmd == "dummy":
        from domtblop_dummy import run
        run(sys.argv[2:])

    elif cmd == "addseq":
        from domtblop_addseq import run
        run(sys.argv[2:])

    else:
        print(__doc__)
        print(f"[domtblop.py] Error: Unknown operation '{cmd}'", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()