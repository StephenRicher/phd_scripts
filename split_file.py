#!/usr/bin/env python3


''' Top-level docstring. Script description. Can be accessed anywhere with 
__doc__ which can be supplied to the argparse description. 
'''

import argparse, sys, logging, gzip, contextlib, binascii, re, stat, os, select

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
        'infiles', nargs = '*', default = ['-'],
        help = 'Specify input files.')
    parser.add_argument(
        '-d', '--gunzip', action = 'store_true',
        help = 'Uncompress GZIP input.')
    parser.add_argument(
        '-c', '--compress', action = 'store_true',
        help = 'GZIP compress ouput.')
    parser.add_argument(
        '-l', '--logfile', nargs = '?', default = None,
        help = 'Specify log file name')
    parser.add_argument(
        '-v', '--verbose', action = 'store_true',
        help = 'Verbose output for debugging.')
    parser.add_argument(
        '-p', '--prefix', default = None,
        type = str,  
        help = 'Specify output file prefix.')    
    parser.add_argument(
        '-n', '--lines', default = 1000,
        type = positive_int,  
        help = 'Specify number of lines per file.')

    args = parser.parse_args()
    
    data_in_stdin = select.select([sys.stdin,],[],[],0.0)[0]
    if not data_in_stdin and args.infiles == ['-']:
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

def main(infiles, logfile, verbose, gunzip, compress, prefix, lines):
    
    ''' Iterate through each infile. '''
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')

    for infile in infiles:
        if not prefix:
            if infile == '-' or named_pipe(infile):
                log.info(
                    f'Unable to set appropriate file prefix - setting to "x".')
                prefix = 'x'
            else:
                prefix = infile
        with gzip_read(infile, gunzip) as i:
            log.info(f'Reading from {infile}')
            counter = 0
            sub = 1
            lines_remaining = True
            while lines_remaining:
                out_name = f'{prefix}{sub}'
                try:
                    if compress:
                        out_obj = gzip.open(f'{out_name}.gz', 'wt')
                    else:
                        out_obj = open(f'{out_name}', 'wt')
                except IOError:
                    log.exception(f'Unable to open {out_name}.')
                    sys.exit()
                with out_obj as w:
                    while counter < lines:
                        counter += 1
                        try:
                            w.write(next(i))
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
        if args['logfile'] == '-':
            (logging.basicConfig(
                stream = sys.stdout,
                format = log_format, 
                level = log_level))
        else:
            (logging.basicConfig(
                filename = args['logfile'],
                format = log_format, 
                level = log_level))
        sys.excepthook = handle_exception
        main(**args)
    finally:
        logging.shutdown()
