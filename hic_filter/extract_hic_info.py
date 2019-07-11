#!/usr/bin/env python3

### Script to extract the HiC read pair information encoding in the SAM optional field 
### by hic_filter.py. Data is output to stdout in long format for analysis.

import sys, argparse, logging

sys.path.insert(0, '/home/stephen/h/phd/scripts2/hic_scripts/hic_filter/')

from sam_class import *
from hic_filter_functions import *

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--sample', help = 'Specify sample name.', type = str, default = None)
parser.add_argument('-f', '--file', nargs = '?', help = 'Specify sam file.', type = argparse.FileType('r'), default = sys.stdin)
parser.add_argument('-o', '--out', nargs = '?', help = 'Specify output file name', type = argparse.FileType('w'),  default = sys.stdout)
parser.add_argument('-l', '--log', nargs = '?', help = 'Specify log file name', type = argparse.FileType('w'),  default = sys.stderr)
args = parser.parse_args()

def handle_exception(exc_type, exc_value, exc_traceback):
    ''' Handle uncaught exceptions (excluding KeyboardInterrupt) 
        through the logging module. '''
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logging.getLogger(__name__).critical("Uncaught exception", exc_info = (exc_type, exc_value, exc_traceback))

sys.excepthook = handle_exception

logging.basicConfig(stream = args.log, format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s', level = logging.DEBUG)
log = logging.getLogger(__name__)

def run():
    try:
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
                        log.exception('Odd number of alignments in file')
                        sys.exit(1)
                    if not is_valid(read1, read2):
                        log.error(f'Invalid format in {read1.qname}.')
                    out.write(f"{args.sample}\t{read1.optional['or:Z']}\t{read1.optional['it:Z']}\t{read1.optional['dt:i']}\t{read1.optional['is:i']}\t{read1.optional['fs:i']}\n")
    finally:
        logging.shutdown()
        
if __name__ == "__main__":
    run()

