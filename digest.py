#!/usr/bin/env python3

import re, argparse, sys, gzip, contextlib, os, stat, select

def restriction_seq(value):
    ''' Custom argparse type restriction enzyme type. '''
    if value.count('^') != 1:
        raise argparse.ArgumentTypeError(f'Error: Restriction site {value} must contain one "^" to indicate cut site.')
    elif re.search('[^ATCG^]', value, re.IGNORECASE)
        raise argparse.ArgumentTypeError(f'Error: Restriction site {value} must only contain "ATCG^".')
    else:
        return value

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', nargs = '?', help = 'Specify input file name.', type = argparse.FileType('r'), default = sys.stdin)
parser.add_argument('-o', '--out', nargs = '?', help = 'Specify output file name', type = argparse.FileType('w'),  default = sys.stdout)
parser.add_argument('-l', '--log', nargs = '?', help = 'Specify log file name', type = argparse.FileType('w'),  default = sys.stderr)
# Create group for required named arguments.
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-r', '--restrictionSite', type = restriction_seq, help = 'Specify restriction cut site (e.g. Mbo1: ^GATC).', required = True)
args = parser.parse_args()

if __name__ == '__main__':
    if not args.file and sys.stdin.isatty():
        parser.print_help()
    
    overhang = args.restrictionSite.index('^')
    site = args.restrictionSite.replace('^','')

    with args.file as f, args.out as out, args.log as log:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
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
                        out.write(f'{ref}\t{start}\t{end}\t{index}\n')    
                    previous_end = end
                out.write(f'{ref}\t{previous_end + 1}\t{len(line)}\t{index + 1}\n')
