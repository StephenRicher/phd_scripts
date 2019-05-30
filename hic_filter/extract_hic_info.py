#!/usr/bin/env python3

### Script to extract the HiC read pair information encoding in the SAM optional field 
### by hic_filter.py. Data is output to stdout in long format for analysis.

import argparse

from sam_class import *
from smart_read import *

parser = argparse.ArgumentParser()
parser.add_argument('file', type = argparse.FileType('r'), help = "Specify sam file.")
args = parser.parse_args()

def main():

    with smart_read(args.file) as sam_file:
        sys.stdout.write('orientation\tinteraction_type\tditag_length\tinsert_size\tfragment_seperation\n')
        for line in sam_file:
            if line.startswith("@"):
                continue
            else:
                try:
                    read1 = sam(line.split())
                    read2 = sam(next(sam_file).split())
                except StopIteration:
                    sys.exit("Error: odd number of alignments in file")
                if not is_valid(read1, read2):
                    raise ValueError(f'Invalid format in {read1.qname}.') 
                sys.stdout.write(f"{read1.optional['or:Z']}\t{read1.optional['it:Z']}\t{read1.optional['dt:i']}\t{read1.optional['is:i']}\t{read1.optional['fs:i']}\n")

if __name__ == "__main__":
    main()

