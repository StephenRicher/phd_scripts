#!/usr/bin/env python3

""" Filter HiC read pairs to remove potential sources of contamination.
    Input should be in SAM/BAM format and have been processed by 
    hictools process.
    """

import sys, argparse, logging

from sam_class import *
from hic_filter_functions import *
from sam_opener import *

def description():
    
    ''' Returns top-level docstring. Useful for providing
        descriptions to sub-parsers after importing.
        '''
        
    return __doc__

def filter(
    infile, output, samtools, sam_out, min_inward, min_outward, max_ditag):

    ''' Iterate through each infile. '''
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')

    if min_inward == min_outward == max_ditag == None:
        log.error('No filter settings defined.')
        sys.exit(1)

    mode = 'wt' if sam_out else 'wb'

    with sam_open(output, mode, samtools = samtools) as out_obj, \
            sam_open(infile, samtools = samtools) as in_obj:
        log.info(f'Writing output to {output}.')

        for line in in_obj:
            if line.startswith("@"):
                out_obj.write(line)
            else:
                try:
                    read1 = sam(line.split())
                    read2 = sam(next(in_obj).split())
                except StopIteration:
                    log.exception('Odd number of alignments in file.')
                    sys.exit(1)
                if not is_valid(read1, read2):
                    log.error(f'Invalid format in {read1.qname}.') 
                elif read1.optional['it:Z'] == "cis":
                    if read1.optional['fs:i'] == 0:
                        continue
                    if max_ditag is not None:
                        if read1.optional['dt:i'] > max_ditag:
                            continue
                    if read1.optional['or:Z'] == 'Inward':
                        if min_inward is not None:
                            if read1.optional['is:i'] < min_inward:
                                continue
                    elif read1.optional['or:Z'] == 'Outward':
                        if min_outward is not None:
                            if read1.optional['is:i'] < min_outward:
                                continue
                out_obj.write(read1.print_sam())
                out_obj.write(read2.print_sam())
