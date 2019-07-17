#!/usr/bin/env python3


''' Top-level docstring. Script description. Can be accessed anywhere with 
__doc__ which can be supplied to the argparse description. 
'''

import argparse, sys, logging, re, stat, os, select
from functools import partial

sys.path.insert(0, '/home/stephen/phd/scripts/python3_modules/')
from gzip_opener import *
from exception_logger import *

def _args():
    
    def positive_int(value):
        
        ''' Custom argparse type to force positive integer. '''
        
        ivalue = int(value)
        if ivalue <= 0:
            raise argparse.ArgumentTypeError(
                f'{value} is not a positive integer.')
        return ivalue

    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'infile', nargs = '?', default = '-',
        help = 'Specify input files.')
    parser.add_argument(
        '-p', '--prefix', default = None,
        type = str,  
        help = 'Specify output file prefix.')    
    parser.add_argument(
        '-n', '--lines', default = 1000,
        type = positive_int,  
        help = 'Specify number of lines per file.')
    parser.add_argument(
        '-l', '--log', nargs = '?', default = None,
        help = 'Specify log file name')
    parser.add_argument(
        '-v', '--verbose', action = 'store_true',
        help = 'Verbose output for debugging.')

    args = parser.parse_args()
    
    data_in_stdin = select.select([sys.stdin,],[],[],0.0)[0]
    if not data_in_stdin and args.infile == '-':
        sys.stderr.write(f'No input provided.\n\n')
        parser.print_help()
        sys.exit(1)
    
    return vars(args)

def named_pipe(path):
    
    """ Check if file is a named pipe. """
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')
    
    if stat.S_ISFIFO(os.stat(path).st_mode):
        pipe = True
    else:
        pipe = False
    return pipe

def main(infile, prefix, lines):
    
    ''' Iterate through each infile. '''
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')

    if not prefix:
        if infile == '-' or named_pipe(infile):
            log.info(
                f'Unable to set appropriate file prefix - setting to "x".')
            prefix = 'x'
        else:
            prefix = in_path
            
    with smart_open(infile, 'r') as in_obj:
        log.info(f'Reading from {infile}')
        counter = 0
        sub = 1
        lines_remaining = True
        while lines_remaining:
            output = f'{prefix}{sub}'
            try:
                out_obj = open(f'{out_path}', 'wt')
            except IOError:
                log.exception(f'Unable to open {out_path}.')
                sys.exit(1)
            with out_obj:
                while counter < lines:
                    counter += 1
                    try:
                        out_obj.write(next(in_obj))
                    except StopIteration:
                        lines_remaining = False
                if lines_remaining:
                    counter = 0
                    sub += 1

if __name__ == '__main__':
    try:
        args = _args()
        fun_name = sys._getframe().f_code.co_name
        log = logging.getLogger(f'{__name__}.{fun_name}')
        
        log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        log_level = logging.DEBUG if args['verbose'] else None
        (logging.basicConfig(
            filename = args['logfile'],
            format = log_format, 
            level = log_level))
        sys.excepthook = handle_exception
        
        # Remove iterator argument and those not needed in main().
        [args.pop(key) for key in ['verbose', 'log']]
        
        main(**args)
        
    finally:
        logging.shutdown()

