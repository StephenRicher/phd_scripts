#!/usr/bin/env python3

### Filter HiC read pairs based on the HiC read pair information encoded in SAM optional field

import argparse, sys

from sam_class import *
from smart_read import *
from hic_filter_functions import *

def positive_int(value):
    """ Custom argparse type for positive integer. """
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(f'{value} is not a positive integer.')
    return ivalue

parser = argparse.ArgumentParser()
parser.add_argument('--min_inward', default = None, type = positive_int, help = 'Specify mininum insert size for inward facing read pairs.')
parser.add_argument('--min_outward', default = None, type = positive_int, help = 'Specify mininum insert size for outward facing read pairs.')
parser.add_argument('--max_ditag', default = None, type = positive_int, help = 'Specify maximum ditag size for read pairs.')
parser.add_argument('file', type = argparse.FileType('r'), help = "Specify sam file.")
args = parser.parse_args()

def main():
    if args.min_inward == args.min_outward == args.max_ditag == None:
        sys.exit('No filter settings defined.')
    
    with smart_read(args.file) as sam_file:
        for line in sam_file:
            if line.startswith("@"):
                sys.stdout.write(line)
            else:
                try:
                    read1 = sam(line.split())
                    read2 = sam(next(sam_file).split())
                except StopIteration:
                    sys.exit("Error: odd number of alignments in file")
                if not is_valid(read1, read2):
                    raise ValueError(f'Invalid format in {read1.qname}.') 
                elif read1.optional['it:Z'] == "cis":
                    if read1.optional['fs:i'] == 0:
                        continue
                    if args.max_ditag is not None:
                        if read1.optional['dt:i'] > args.max_ditag:
                            continue
                    if read1.optional['or:Z'] == "Inward":
                        if args.min_inward is not None:
                            if read1.optional['is:i'] < args.min_inward:
                                continue
                    elif read1.optional['or:Z'] == "Outward":
                        if args.min_outward is not None:
                            if read1.optional['is:i'] < args.min_outward:
                                continue
                sys.stdout.write(read1.print_sam())
                sys.stdout.write(read2.print_sam())

if __name__ == "__main__":
    main()

