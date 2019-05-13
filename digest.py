#!/usr/bin/env python3

import re, argparse, sys, gzip, contextlib, os, stat, select

from smart_open import smart_open

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help = "Specify input file name.", nargs = "?")
parser.add_argument("-o", "--out", help = "Specify output file name", nargs = "?")
parser.add_argument("-r", "--restrictionSite", help = "Specify restriction cut site (e.g. Mbo1: ^GATC).")
args = parser.parse_args()

def eprint(*args, **kwargs):
    """ Print to stderr and exit """
    print(*args, file=sys.stderr, **kwargs)
    sys.exit(1)

if __name__ == '__main__':
    if not args.file and sys.stdin.isatty():
        parser.print_help()

    if args.restrictionSite.count("^") != 1:
        eprint("Error: Restriction site {} must contain one '^' to indicate cut site.".format(args.restrictionSite))
    
    overhang = args.restrictionSite.index('^')
    site = args.restrictionSite.replace('^','')
    
    if re.search("[^ATCG]", site, re.IGNORECASE):
        eprint("Error: Restriction sequence {} must only contain ATCG.".format(site))
    
    with smart_open(args.file, mode = "read") as f, smart_open(args.out, mode = "write") as out:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith(">"):
                ref = line.rsplit()[0][1:]
            else:
                matches = re.finditer(site, line)
                index = 0
                previous_end = 0
                for match in matches:
                    # Skip if restriction sequence at start of reference.
                    if match.start() != 0:
                        index += 1
                        start = 1 if index == 1 else previous_end + 1
                        end = match.start() + overhang
                        out.write("{}\t{}\t{}\t{}\n".format(ref, start, end, index))    
                    previous_end = end
                out.write("{}\t{}\t{}\t{}\n".format(ref, previous_end + 1, len(line), index + 1))
