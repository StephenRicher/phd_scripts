#!/usr/bin/env python3

""" Extract HiC read pair information encodiing by hictools process.
    Useful for assessing contamination and determining parameters for
    hictools filter.
"""

import sys, argparse, logging

from sam_class import *
from hic_filter_functions import *
from sam_opener import *
from gzip_opener import *

def description():
    
    ''' Returns top-level docstring. Useful for providing
        descriptions to sub-parsers after importing.
        '''
        
    return __doc__

def extract(
    infile, output, samtools, sample, write_gzip):
    if not sample:
        sample = infile
        
    with smart_open(output, 'wt', write_gzip) as out_obj, \
            sam_open(infile, samtools = samtools) as in_obj:
        log.info(f'Writing output to {output}.')
        out_obj.write(
            'sample\torientation\tinteraction_type\tditag_length\t'
            'insert_size\tfragment_seperation\n')
        for i, line in enumerate(in_obj):
            if line.startswith('@'):
                continue
            else:
                try:
                    read1 = sam(line.split())
                    read2 = sam(next(in_obj).split())
                except StopIteration:
                    log.exception('Odd number of alignments in file')
                    sys.exit(1)
                if not is_valid(read1, read2):
                    log.error(f'Invalid format in {read1.qname}.')
                out_obj.write(
                    f'{sample}\t{read1.optional["or:Z"]}\t'
                    f'{read1.optional["it:Z"]}\t{read1.optional["dt:i"]}\t'
                    f'{read1.optional["is:i"]}\t{read1.optional["fs:i"]}\n')

