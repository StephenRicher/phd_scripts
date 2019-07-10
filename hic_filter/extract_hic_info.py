#!/usr/bin/env python3

### Script to extract the HiC read pair information encoding in the SAM optional field 
### by hic_filter.py. Data is output to stdout in long format for analysis.

import argparse

sys.path.insert(0, '/home/stephen/h/phd/scripts2/hic_scripts/hic_filter/')

from sam_class import *
from hic_filter_functions import *

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--sample', help = 'Specify sample name.', type = str, default = None)
parser.add_argument('-f', '--file', nargs = '?', help = 'Specify sam file.', type = argparse.FileType('r'), default = sys.stdin)
parser.add_argument('-o', '--out', nargs = '?', help = 'Specify output file name', type = argparse.FileType('w'),  default = sys.stdout)
args = parser.parse_args()

def main():
    
    if not args.sample:
        if args.file is not sys.stdin:
            args.sample = args.file
        else:
            args.sample = 'sample'
        
    with args.file as sam_file, args.out as out:
        out.write('sample\torientation\tinteraction_type\tditag_length\tinsert_size\tfragment_seperation\n')
        for i, line in enumerate(sam_file):
            if line.startswith('@'):
                continue
            else:
                try:
                    read1 = sam(line.split())
                    read2 = sam(next(sam_file).split())
                except StopIteration:
                    sys.exit('Error: odd number of alignments in file')
                if not is_valid(read1, read2):
                    raise ValueError(f'Invalid format in {read1.qname}.') 
                out.write(f"{args.sample}\t{read1.optional['or:Z']}\t{read1.optional['it:Z']}\t{read1.optional['dt:i']}\t{read1.optional['is:i']}\t{read1.optional['fs:i']}\n")
                
if __name__ == "__main__":
    main()

