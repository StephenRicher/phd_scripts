#!/usr/bin/env python3

import re, argparse, sys, gzip, contextlib, os, stat, select

from smart_open import smart_open

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help = "Specify input file name.", nargs = "?")
args = parser.parse_args()

def eprint(*args, **kwargs):
    """ Print to stderr and exit """
    print(*args, file=sys.stderr, **kwargs)
    sys.exit(1)

if __name__ == '__main__':
    if not args.file and sys.stdin.isatty():
        parser.print_help()

    with smart_open(args.file, mode = "read") as f:
        
        for index, line in enumerate(f):
            line = line.rstrip('\n')
            in_header = True if line.startswith(">") else False
            header_above = None
            if not line:
                eprint("Error: Line {} is a blank line.".format(index))
            elif not in_header:
                invalid_characters = re.search("[^ATCGURYKMSWBDHVN-]", line, re.IGNORECASE)
                if invalid_characters:
                    eprint("Error: Invalid sequence character '{}' at line {}.".format(invalid_characters.group(0), index))
                elif in_header == header_above:
                    eprint("Error: FASTA sequence contains line wrapping.")
            # If we are in header and the previous line was also a header...
            elif in_header == header_above:
                eprint("Error: Header missing sequence at line {}.".format(index))
            header_above = in_header
        if in_header:
            eprint("Error: Header missing sequence on last last.")
        
