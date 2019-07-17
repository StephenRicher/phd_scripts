#!/usr/bin/env python3

### Filter HiC read pairs based on the HiC read pair information encoded in SAM optional field

import argparse, sys, logging

sys.path.insert(0, '/home/stephen/h/phd/scripts2/hic_scripts/hic_filter/')

from sam_class import *
from hic_filter_functions import *
from gzip_opener import *

def positive_int(value):
    ''' Custom argparse type for positive integer. '''
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(f'{value} is not a positive integer.')
    return ivalue

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', nargs = '?', help = 'Specify sam file.', type = argparse.FileType('r'), default = sys.stdin)
parser.add_argument('-o', '--out', nargs = '?', help = 'Specify output file name', type = argparse.FileType('w'),  default = sys.stdout)
parser.add_argument('-l', '--log', nargs = '?', help = 'Specify log file name', type = argparse.FileType('w'),  default = sys.stderr)
parser.add_argument('--min_inward', default = None, type = positive_int, help = 'Specify mininum insert size for inward facing read pairs.')
parser.add_argument('--min_outward', default = None, type = positive_int, help = 'Specify mininum insert size for outward facing read pairs.')
parser.add_argument('--max_ditag', default = None, type = positive_int, help = 'Specify maximum ditag size for read pairs.')
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
        if args.min_inward == args.min_outward == args.max_ditag == None:
            log.error('No filter settings defined.')
            sys.exit(1)
    
        with args.file as sam_file, args.out as out:
            for line in sam_file:
                if line.startswith("@"):
                    out.write(line)
                else:
                    try:
                        read1 = sam(line.split())
                        read2 = sam(next(sam_file).split())
                    except StopIteration:
                        log.exception('Error: odd number of alignments in file.')
                        sys.exit(1)
                    if not is_valid(read1, read2):
                        log.error(f'Invalid format in {read1.qname}.') 
                    elif read1.optional['it:Z'] == "cis":
                        if read1.optional['fs:i'] == 0:
                            continue
                        if args.max_ditag is not None:
                            if read1.optional['dt:i'] > args.max_ditag:
                                continue
                        if read1.optional['or:Z'] == 'Inward':
                            if args.min_inward is not None:
                                if read1.optional['is:i'] < args.min_inward:
                                    continue
                        elif read1.optional['or:Z'] == 'Outward':
                            if args.min_outward is not None:
                                if read1.optional['is:i'] < args.min_outward:
                                    continue
                    out.write(read1.print_sam())
                    out.write(read2.print_sam())
    finally:
        logging.shutdown()

if __name__ == '__main__':
    run()

