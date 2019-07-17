#!/usr/bin/env python3

'''
Randomly sample lines from a file. Blocksize can be specified to 
sample lines in blocks. This is useful for sampling FASTQ files where
each FASTQ record consists of 4 line blocks.
'''

import argparse, sys, logging, gzip, contextlib, binascii, re, random, itertools, select, multiprocessing
from functools import partial
from itertools import islice, chain

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
        help = 'Input file name.')
    parser.add_argument(
        '-o', '--output', nargs = '?', default = '-', 
        help = 'Output file name.')
    parser.add_argument(
        '-n', '--nsamples', default = 1,
        type = positive_int,  
        help = 'Specify sample size.')
    parser.add_argument(
        '-b', '--blocksize', default = 1,
        type = positive_int,  
        help = 'Specify number lines per block.')
    parser.add_argument(
        '-l', '--log', nargs = '?', default = None,
        help = 'Specify log file name')
    parser.add_argument(
        '-v', '--verbose', action = 'store_true',
        help = 'Verbose output for debugging.')

    args = parser.parse_args()
    
    data_in_stdin = select.select([sys.stdin,],[],[],0.0)[0]
    if not data_in_stdin and args.infile == '-':
        sys.stderr.write(f'Error: No input provided.\n\n')
        parser.print_help()
        sys.exit(1)
    
    return vars(args)

def reservoir_sample(in_obj, nsample, blocksize):
    """ Perform reservoir sampling on lines in file object.
    
    Return 'nsample' random blocks where each block is 'blocksize' lines.
    https://stackoverflow.com/questions/2612648/reservoir-sampling
    """
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')
    
    reservoir = []
    t = -1
    nlines = 0
    while True:
        block = [x.strip('\n') 
                for x 
                in list(itertools.islice(in_obj, 0, blocksize))]
        t += 1
        nlines += len(block)
        if not block:
            break
        elif t < nsample:
            reservoir.append(block)
        else:
            m = random.randint(0, t)
            if m < nsample:
                reservoir[m] = block
    
    if nlines % blocksize != 0:
        msg = f'Total lines ({nlines}) not divisible by blocksize ({blocksize}).'
        log.error(msg)
        
    # Flatten list of lists
    reservoir = list(itertools.chain(*reservoir))
    return reservoir
    
def main(infile, output, nsamples, blocksize):
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')

    with smart_open(output, 'w') as out_obj, \
            smart_open(infile, 'r') as in_obj:
        log.info(f'Writing output to {output}')
        sample = reservoir_sample(in_obj, nsamples, blocksize)
        for line in sample:
            out_obj.write(f'{line}\n')

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
