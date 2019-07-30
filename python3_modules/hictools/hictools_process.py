#!/usr/bin/env python3

""" Process named-sorted SAM/BAM alignment files to identify fragment 
    mapping, insert size, ditag size and relative orientation of pairs.
"""

import argparse, collections, bisect, time, re, sys, math, logging

from sam_class import *
from hic_filter_functions import *
from sam_opener import *
from gzip_opener import *

def description():
    
    ''' Returns top-level docstring. Useful for providing
        descriptions to sub-parsers after importing.
        '''
        
    return __doc__

def process(
    infile, output, read_gzip, samtools, sam_out, digest):
    
    mode = 'wt' if sam_out else 'wb'
        
    with sam_open(output, mode, samtools = samtools) as out_obj, \
            sam_open(infile, samtools = samtools) as in_obj, \
            smart_open(digest, 'rt', read_gzip) as digest:
        d = process_digest(digest)
        for line in in_obj:
            if line.startswith("@"):
                out_obj.write(line)
                continue
            else:
                try:
                    read1 = sam(line.split())
                    read2 = sam(next(in_obj).split())
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
                out_obj.write(read1.print_sam())
                out_obj.write(read2.print_sam())


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

