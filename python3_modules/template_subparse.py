#!/usr/bin/env python3


''' Top-level docstring. Script description. Can be accessed anywhere with 
__doc__ which can be supplied to the argparse description This template. 
'''

import argparse, sys, logging, select

from exception_logger import *
from gzip_opener import *

def main():

    epilog = 'Steven Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'

    formatter_class = argparse.ArgumentDefaultsHelpFormatter

    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = formatter_class,
        epilog = epilog)
        
    # Set parent parser for options common to all sub commands
    base_parser = argparse.ArgumentParser(
        formatter_class = formatter_class,
        add_help = False)
    base_parser.add_argument(
        'infile', nargs = '?', default = '-',
        help = 'Input file name.')
    base_parser.add_argument(
        '-o', '--output', nargs = '?', default = '-', 
        help = 'Output file name.')
    base_parser.add_argument(
        '-z', '--gzip', 
        action = 'store_true', dest = 'write_gzip',
        help = 'Compress output using gzip')
    base_parser.add_argument(
        '-u', '--gunzip', 
        action = 'store_true', dest = 'read_gzip',
        help = 'Read gzip compressed input.')
    base_parser.add_argument(
        '-l', '--log', nargs = '?',
        help = 'Specify log file name.')
    base_parser.add_argument(
        '-v', '--verbose', 
        action = 'store_true',
        help = 'number')

    # Define subparser
    subparsers = parser.add_subparsers(
        title = 'required commands',
        description = '',
        dest = 'command',
        metavar = 'Commands',
        help = 'Description:')
    
    # Store each command name and parser in dictionary
    commands = {}
    
    # Cat sub-parser
    cat_command = 'cat'
    cat_parser = subparsers.add_parser(cat_command,
        description = 'cat command help description',
        help = 'Write input data to output.', 
        parents = [base_parser],
        formatter_class = formatter_class,
        epilog = epilog)
    cat_parser.add_argument(
        '-p', '--parameter', default = 50,
        type = positive_int,
        help = 'number')
    cat_parser.set_defaults(function = cat)
    commands[cat_command] = cat_parser
    
    # Rev sub-parser
    rev_command = 'rev'
    rev_parser = subparsers.add_parser(rev_command, 
        description = 'rev command help description',
        help = 'Reverse input data and write to output.', 
        parents = [base_parser],
        formatter_class = formatter_class,
        epilog = epilog)
    rev_parser.add_argument(
        '-n', '--number', default = 20,
        type = int,
        help ='number')
    rev_parser.set_defaults(function = rev)
    commands[rev_command] = rev_parser
    
    args = parser.parse_args()
    
    try:
        func = args.function
    except AttributeError:
        parser.print_help()
        sys.exit(1)
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    log_level = logging.DEBUG if args.verbose else None
    (logging.basicConfig(
        filename = args.log,
        format = log_format, 
        level = log_level))
    sys.excepthook = handle_exception

    data_in_stdin = select.select([sys.stdin,],[],[],0.0)[0]
    if not data_in_stdin and args.infile == '-':
        log.error(f'No input provided.\n')
        # Call specific parser help for evoked command
        commands[args.command].print_help()
        sys.exit(1)

    args_dict = vars(args)
    [args_dict.pop(key) for key in ['command', 'function', 'verbose', 'log']]
    return func(**vars(args))

def positive_int(value):
        
    ''' Custom argparse type to force positive integer. '''
        
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(
            f'{value} is not a positive integer.')
    return ivalue

def cat(infile, output, read_gzip, write_gzip, parameter):
    
    ''' Iterate through each infile. '''

    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')
    
    with smart_open(output, mode = 'wt', gz = write_gzip) as out_obj, \
            smart_open(infile, mode = 'rt', gz = read_gzip) as in_obj:
        for line in in_obj:
            out_obj.write(line)

def rev(infile, output, read_gzip, write_gzip, number):
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')
    
    with smart_open(output, mode = 'wt', gz = write_gzip) as out_obj, \
            smart_open(infile, mode = 'rt', gz = read_gzip) as in_obj:
        for line in in_obj:
            out_obj.write(line[::-1])

if __name__ == '__main__':
    try:
        sys.exit(main())
    finally:
        logging.shutdown()


