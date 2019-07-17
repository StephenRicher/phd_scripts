#!/usr/bin/env python3


''' Truncate proximity ligated restriction digest fragments at the 
    ligation junction.
'''

import os, argparse, sys, logging, gzip, contextlib, binascii, re, select, time

from gzip_opener import *
from exception_logger import *

def description():
    
    ''' Returns top-level docstring. Useful for providing
        descriptions to sub-parsers after importing.
        '''
        
    return __doc__

def process_restriction(restriction):
    restriction = restriction.upper()
    restriction_seq = restriction.replace('^','')
    cut_site1 = restriction.index('^')
    ligation2 = restriction_seq[cut_site1:]
    cut_site2 = len(restriction) - cut_site1 - 1
    ligation1 = restriction_seq[0:cut_site2]
    ligation_seq = ligation1 + ligation2
    cut_index = len(ligation1)
    return ligation_seq, restriction_seq

def truncate(infile, output, read_gzip, write_gzip, sample, summary, restriction):
    
    ''' Run main loop. '''

    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')

    ligation_seq, restriction_seq = process_restriction(restriction)
    truncated = 0
    truncated_length = 0

    if not sample:
        sample = infile

    with smart_open(output, 'wt', write_gzip) as out_obj, \
            smart_open(infile, 'rt', read_gzip) as in_obj, \
            open(summary, 'a') as summary_obj:
            
        is_truncated = False
        for index, line in enumerate(in_obj):
            line=line.rstrip('\n')
            # Sequence line
            if index % 4 == 1:
                line = line.upper()
                if ligation_seq in line:
                    line = line[0: line.index(ligation_seq)] + restriction_seq
                    is_truncated = True
                seq_length = len(line)
            # Quality line
            elif index % 4 == 3:
                line = line[0:seq_length]
                if is_truncated:
                    truncated += 1
                    truncated_length += seq_length
                    is_truncated = False
            out_obj.write(f'{line}\n')
        try:
            mean_truncated_length = truncated_length/truncated
        except ZeroDivisionError:
            mean_truncated_length = 'na'
        summary_obj.write(
            f'{sample}\t{truncated}\t{ mean_truncated_length}\n')
