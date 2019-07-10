#!/usr/bin/env python3

import argparse, sys, re

def restriction_seq(value):
    ''' Custom argparse type restriction enzyme type. '''
    if value.count('^') != 1:
        raise argparse.ArgumentTypeError(f'Error: Restriction site {value} must contain one "^" to indicate cut site.')
    elif re.search('[^ATCG^]', value, re.IGNORECASE)
        raise argparse.ArgumentTypeError(f'Error: Restriction site {value} must only contain "ATCG^".')
    else:
        return value

parser = argparse.ArgumentParser(description = 'Script to truncate proximity ligated restriction fragments in FASTQ files.')
parser.add_argument('-f', '--file', nargs = '?', help = 'Specify input file name.', type = argparse.FileType('r'), default = sys.stdin)
parser.add_argument('-o', '--out', nargs = '?', help = 'Specify output file name', type = argparse.FileType('w'),  default = sys.stdout)
parser.add_argument('-l', '--log', nargs = '?', help = 'Specify log file name', type = argparse.FileType('w'),  default = sys.stderr)
parser.add_argument('-s', '--sample', help = 'Specify sample name.', type = str, default = None)
# Create group for required named arguments.
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-r', '--restrictionSite', type = restriction_seq, help = 'Specify restriction cut site (e.g. Mbo1: ^GATC).', required = True)
args = parser.parse_args()

def process_restriction(restriction):
    restriction_seq = restriction.replace('^','')
    cut_site1 = restriction.index('^')
    ligation2 = restriction_seq[cut_site1:]
    cut_site2 = len(restriction) - cut_site1 - 1
    ligation1 = restriction_seq[0:cut_site2]
    ligation_seq = ligation1 + ligation2
    cut_index = len(ligation1)
    return ligation_seq, restriction_seq

if __name__ == '__main__':
    if not args.file and sys.stdin.isatty():
        parser.print_help()
 
    if not args.sample:
        if args.file is not sys.stdin:
            args.sample = args.file
        else:
            args.sample = 'sample'
 
    ligation_seq, restriction_seq = process_restriction(args.restrictionSite)
    truncated = 0
    truncated_length = 0
    
    with args.file as f, args.out as out, args.log as log:
        for index, line in enumerate(f):
            line=line.rstrip('\n')
            # Sequence line
            if index % 4 == 1 and ligation_seq in line:
                line = line[0: line.index(ligation_seq)] + restriction_seq
                is_truncated = True
            seq_length = len(line)
            # Quality line
            elif index % 4 == 3:
                line = line[0:seq_length]
                if truncated:
                    truncated += 1
                    truncated_length += seq_length
                    is_truncated = False
            out.write(f'{line}\n')
        log.write(f'{args.sample}\t{truncated}\t{truncated_length/truncated}\n')
            
