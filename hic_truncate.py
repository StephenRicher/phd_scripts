#!/usr/bin/env python3

import re, argparse, sys, gzip, contextlib, os, stat, select

from smart_open import smart_open

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help = "Specify input file name.", nargs = "?")
parser.add_argument("-r", "--restriction", help = "Specific restriction enzyme with cut site (i.e. Mbol: ^GATC)")
args = parser.parse_args()

@contextlib.contextmanager
def smart_read (filename):
    """ Custom context file opener which can read from from stdin or input file."""
    
    if filename and filename != '-':
        try:
            fh = open(filename, 'rt')
        except IOError as e:
            sys.exit(str(e))
    else:
        if not sys.stdin.isatty():
            fh = sys.stdin
        else:
            sys.exit("Nothing in stdin and no infile provided.")
    try:
        yield fh
    finally:
        if fh not in [sys.stdin]:
            fh.close()

def process_restriction(restriction):
    restriction_seq = restriction.replace('^','')
    cut_site1 = restriction.index('^')
    ligation2 = restriction_seq[cut_site1:]
    cut_site2 = len(restriction) - cut_site1 - 1
    ligation1 = restriction_seq[0:cut_site2]
    ligation_seq = ligation1 + ligation2
    cut_index = len(ligation1)
    return ligation_seq, restriction_seq

def eprint(*args, **kwargs):
    """ Print to stderr """
    print(*args, file = sys.stderr, **kwargs)

if __name__ == '__main__':
    if not args.file and sys.stdin.isatty():
        parser.print_help()
 
    ligation_seq, restriction_seq = process_restriction(args.restriction)
    total_truncated = 0
    trunc_len_dist = {}
    
    with smart_read(args.file) as f:
        seq_length = False
        for index, line in enumerate(f):
            line=line.rstrip('\n')
            # Sequence line
            if index % 4 == 1 and ligation_seq in line:
                line = line[0: line.index(ligation_seq)] + restriction_seq
                seq_length = len(line)
            # Quality line
            elif seq_length and index % 4 == 3:
                line = line[0:seq_length]
                if seq_length in trunc_len_dist:
                    trunc_len_dist[seq_length] += 1
                else:
                    trunc_len_dist[seq_length] = 1
                total_truncated += 1
                seq_length = False
            print(line)
    for length, number in trunc_len_dist.items():
        eprint(args.file, length, number)
            
