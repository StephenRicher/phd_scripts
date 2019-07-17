#!/usr/bin/env python3


''' Top-level docstring. Script description. Can be accessed anywhere with 
__doc__ which can be supplied to the argparse description This template. 
'''

import argparse, sys, logging, gzip, contextlib, binascii, re, select, multiprocessing
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
        'infiles', nargs = '*', default = ['-'],
        help = 'Specify input files.')
    parser.add_argument(
        '-o', '--outfiles', nargs = '*', default = ['-'], 
        help = 'Specify output files.')
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
        '-@', '--threads', nargs = '?', default = 1,
        type = positive_int,
        help = 'Number of threads to use in parallel.')
    parser.add_argument(
        '-v', '--verbose', action = 'store_true',
        help = 'Verbose output for debugging.')
    # Create group for required named arguments.
    requiredNamed = parser.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument(
        '-x', '--parameter', required = True,
        type = positive_int,  
        help = 'Specify positive integer.')

    args = parser.parse_args()
    
    data_in_stdin = select.select([sys.stdin,],[],[],0.0)[0]
    if not data_in_stdin and args.infiles == ['-']:
        sys.stderr.write(f'Error: No input provided.\n\n')
        parser.print_help()
        sys.exit(1)
    
    return vars(args)

def main(in_path, out_path, parameter, gunzip, compress):
    
    ''' Iterate through each infile. '''

    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')
    
    with gz_write(out_path, compress) as out_obj, \
            gz_read(in_path, gunzip) as in_obj:
        for line in in_obj:
            out_obj.write(line)


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
        
        # Initiate multiprocessing pool with set number of processes.
        pool = multiprocessing.Pool(processes = args['threads'])
        
        if len(args['infiles']) != len(args['outfiles']): 
            if args['outfiles'] != ['-']:
                log.error('Number of input files does not ' 
                          'match number of output files.')
                sys.exit(1)
            else:
                log.warning('Processing multiple input files to stdout.') 
                args['outfiles'] = args['outfiles'] * len(args['infiles'])

        iter_zip = list(zip(args['infiles'], args['outfiles']))

        # Remove iterator argument and those not needed in main().
        [args.pop(key) for key in 
            ['infiles', 'outfiles', 'verbose', 'logfile', 'threads']]
        
        pool_outputs = pool.starmap(partial(main, **args), iter_zip)
        
    finally:
        pool.close()
        pool.join()
        logging.shutdown()

