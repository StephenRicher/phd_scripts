#!/usr/bin/env python3

import sys, argparse, logging

sys.path.insert(0, '/home/stephen/h/phd/scripts2/hic_scripts/hic_filter/')

from sam_class import *
from hic_filter_functions import *

parser = argparse.ArgumentParser()
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
        with args.file as sam_file, args.out as out:
            a = 0
            for i, line in enumerate(sam_file):
                if line.startswith('@'):
                    continue
                else:
                    try:
                        a += 1
                        read1 = sam(line.split())
                        read2 = sam(next(sam_file).split())
                    except StopIteration:
                        log.exception('Odd number of alignments in file')
                        sys.exit(1)
                    if not is_valid(read1, read2):
                        log.error(f'Invalid format in {read1.qname}.')
                    if a == 1:
                        min_frag = read1.optional['fn:i']
                    str1 = 1 if read1.is_reverse else 0
                    str2 = 1 if read2.is_reverse else 0
                    chr1 = read1.rname
                    chr2 = read2.rname
                    pos1 = read1.left_pos
                    pos2 = read2.left_pos
                    frag1 = read1.optional['fn:i'] - min_frag + 1
                    frag2 = read2.optional['fn:i'] - min_frag + 1
                    out.write(f"{str1}\t{chr1}\t{pos1}\t{frag1}\t{str2}\t{chr2}\t{pos2}\t{frag2}\n")
    finally:
        logging.shutdown()
        
if __name__ == "__main__":
    run()
