#!/usr/bin/env python3

import re, argparse, sys, gzip, contextlib, os, stat, select

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', nargs = '?', help = 'Specify input file name.', type = argparse.FileType('r'), default = sys.stdin)
parser.add_argument('-l', '--log', nargs = '?', help = 'Specify log file name', type = argparse.FileType('w'),  default = sys.stderr)
args = parser.parse_args()

if __name__ == '__main__':
    if not args.file and sys.stdin.isatty():
        parser.print_help()

    invalid = False

    with args.file as f, args.log as log::
        for index, line in enumerate(f):
            line = line.rstrip('\n')
            in_header = True if line.startswith('>') else False
            header_above = None
            if not line:
                log.write(f'Error: Line {index} is a blank line.')
                invalid = True
            elif not in_header:
                invalid_characters = re.search('[^ATCGURYKMSWBDHVN-]', line, re.IGNORECASE)
                if invalid_characters:
                    log.write(f'Error: Invalid sequence character "{invalid_characters.group(0)}" at line {index}.')
                    invalid = True
                elif in_header == header_above:
                    log.write('Error: FASTA sequence contains line wrapping.')
                    invalid = True
            # If we are in header and the previous line was also a header...
            elif in_header == header_above:
                log.write(f'Error: Header missing sequence at line {index}.')
                invalid = True
            header_above = in_header
        if in_header:
            log.write('Error: Header missing sequence on last last.')
            invalid = True
        if invalid:
            sys.exit(1)
        
