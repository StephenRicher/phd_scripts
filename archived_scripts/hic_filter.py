#!/usr/bin/env python3

import argparse, collections, bisect, time, re, sys, math, logging

sys.path.insert(0, '/home/stephen/h/phd/scripts2/hic_scripts/hic_filter/')

from sam_class import *
from hic_filter_functions import *

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', nargs = '?', help = 'Specify sam file.', type = argparse.FileType('r'), default = sys.stdin)
parser.add_argument('-o', '--out', nargs = '?', help = 'Specify output file name', type = argparse.FileType('w'),  default = sys.stdout)
parser.add_argument('-l', '--log', nargs = '?', help = 'Specify log file name', type = argparse.FileType('w'),  default = sys.stderr)
# Create group for required named arguments.
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-d', '--digest', help = 'Specify digest file.', type = argparse.FileType('r'), required = True)
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
        with args.file as sam_file, args.out as out, args.digest as digest:
            d = process_digest(digest)
            for line in sam_file:
                if line.startswith("@"):
                    out.write(line)
                    continue
                else:
                    try:
                        read1 = sam(line.split())
                        read2 = sam(next(sam_file).split())
                    except StopIteration:
                        log.exception("Odd number of alignments in file")
                    orientation, ditag_length, insert_size, interaction, fragment_seperation, read1_fragment, read2_fragment = run_filter(read1, read2, d)
                    read1.optional['or:Z'] = orientation
                    read2.optional['or:Z'] = orientation
                    read1.optional['it:Z'] = interaction
                    read2.optional['it:Z'] = interaction
                    read1.optional['dt:i'] = ditag_length
                    read2.optional['dt:i'] = ditag_length
                    read1.optional['is:i'] = insert_size
                    read2.optional['is:i'] = insert_size
                    read1.optional['fs:i'] = fragment_seperation
                    read2.optional['fs:i'] = fragment_seperation
                    read1.optional['fn:i'] = read1_fragment
                    read2.optional['fn:i'] = read2_fragment
                    out.write(read1.print_sam())
                    out.write(read2.print_sam())
    finally:
        logging.shutdown()
                
def pysam_test():
    import pysam
    samfile = pysam.AlignmentFile("/media/stephen/Data/hic_analysis/data/test.sam", "rb")
    mysamfile = open("/media/stephen/Data/hic_analysis/data/test2.sam")
    for read_pysam, read_mysam in zip(samfile.fetch(), mysamfile):
        myread = sam(read_mysam.split())
        assert read_pysam.get_reference_positions()[-1] + 1 == myread.right_pos
        assert read_pysam.get_reference_positions()[0] + 1 == myread.left_pos
        assert read_pysam.reference_length == myread.reference_length
        assert read_pysam.is_paired == myread.is_paired
        assert read_pysam.is_reverse == myread.is_reverse
        if myread.is_reverse:
            assert read_pysam.get_reference_positions()[-1] + 1 == myread.five_prime_pos
            assert read_pysam.get_reference_positions()[0] + 1 == myread.three_prime_pos
        else:
            assert read_pysam.get_reference_positions()[-1] + 1 == myread.three_prime_pos
            assert read_pysam.get_reference_positions()[0] + 1 == myread.five_prime_pos
    mysamfile.close()
            
if __name__ == "__main__":
    run()
    #pysam_test()
    

