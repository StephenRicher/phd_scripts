#!/usr/bin/env python3


''' Top-level docstring. Script description. Can be accessed anywhere with 
__doc__ which can be supplied to the argparse description This template. 
'''

import os, argparse, sys, logging, gzip, contextlib, binascii, re, select, time

from functools import partial
from gzip_opener import *
from exception_logger import *

def _args():
    
    def restriction_seq(value):
        
        ''' Custom argument type for restriction enzyme argument. '''
        
        if value.count('^') != 1:
            raise argparse.ArgumentTypeError(
                f'Restriction site {value} must contain one "^" at cut site.')
        elif re.search('[^ATCG^]', value, re.IGNORECASE):
            raise argparse.ArgumentTypeError(
                f'Restriction site {value} must only contain "ATCG^".')
        else:
            return value.upper()
    
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'infile', nargs = '?', default = '-',
        help = 'Input file name.')
    parser.add_argument(
        '-o', '--output', nargs = '?', default = '-', 
        help = 'Output file name.')
    parser.add_argument(
        '-s', '--summary', nargs = '?', 
        default = f'truncation_summary_{time.strftime("%Y%m%d-%H%M%S")}.txt', 
        help = 'Truncation summary file name.')
    parser.add_argument(
        '-n', '--sample', default = None,
        help = 'Sample name for truncation summary file.')
    parser.add_argument(
        '-l', '--log', nargs = '?', default = None,
        help = 'Specify log file name')
    parser.add_argument(
        '-v', '--verbose', action = 'store_true',
        help = 'Verbose output for debugging.')
    # Create group for required named arguments.
    requiredNamed = parser.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument(
        '-r', '--restriction', required = True,
        type = restriction_seq, 
        help = 'Restriction cut sequence with "^" to indicate cut site. '
               'e.g. Mbol = ^GATC')

    args = parser.parse_args()
    
    data_in_stdin = select.select([sys.stdin,],[],[],0.0)[0]
    if not data_in_stdin and args.infile == '-':
        sys.stderr.write(f'Error: No input provided.\n\n')
        parser.print_help()
        sys.exit(1)
    
    return vars(args)

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

def main(infile, output, sample, summary, restriction):
    
    ''' Run main loop. '''

    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')

    ligation_seq, restriction_seq = process_restriction(restriction)
    truncated = 0
    truncated_length = 0

    if not sample:
        sample = infile

    with smart_open(infile, 'r') as in_obj, \
            open(summary, 'a') as summary_obj, \
            smart_open(output, 'w') as out_obj:
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

if __name__ == '__main__':
    try:
        args = _args()
        fun_name = sys._getframe().f_code.co_name
        log = logging.getLogger(f'{__name__}.{fun_name}')
        
        log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        log_level = logging.DEBUG if args['verbose'] else None
        (logging.basicConfig(
            filename = args['log'],
            format = log_format, 
            level = log_level))
        sys.excepthook = handle_exception

        # Remove arguments not needed in main().
        [args.pop(key) for key in ['verbose', 'log']]
        
        main(**args)

    finally:
        logging.shutdown()
