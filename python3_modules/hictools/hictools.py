#!/usr/bin/env python3


''' HiCTools is a set of tools for analysing HiC data. '''

import sys, argparse, logging, select, re, time

from exception_logger import *
from gzip_opener import *
import  hictools_digest, hictools_truncate, \
        hictools_filter, hictools_extract, \
        hictools_process

def main():

    epilog = 'Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'

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
        '-o', '--output', nargs = '?', default = '-', 
        help = 'Output file.')
    base_parser.add_argument(
        '-l', '--log', nargs = '?',
        help = 'Log output file.')
    base_parser.add_argument(
        '-v', '--verbose', 
        action = 'store_true',
        help = 'Verbose logging for debugging.')
    
    # Set parent parser for options only non-SAM processing
    txt_parser = argparse.ArgumentParser(
        formatter_class = formatter_class,
        add_help = False)   
    txt_parser.add_argument(
        '-z', '--gzip', 
        action = 'store_true', dest = 'write_gzip',
        help = 'Compress output using gzip')
    txt_parser.add_argument(
        '-u', '--gunzip', 
        action = 'store_true', dest = 'read_gzip',
        help = 'Read gzip compressed input.')
        
    # Set common parent parser options for SAM processing
    sam_parser = argparse.ArgumentParser(
        formatter_class = formatter_class,
        add_help = False)   
    sam_parser.add_argument(
        'infile', nargs = '?', default = '-',
        help = 'Input file in SAM/BAM format.')
    sam_parser.add_argument(
        '--samtools', default = 'samtools',
        help = 'Set path to samtools installation.')

    # Define subparser
    subparsers = parser.add_subparsers(
        title = 'required commands',
        description = '',
        dest = 'command',
        metavar = 'Commands',
        help = 'Description:')
    
    # Dictionary to store each command name and command sub-parser
    commands = {}
    
    # Digest sub-parser
    digest_command = 'digest'
    digest_parser = subparsers.add_parser(digest_command,
        description = hictools_digest.description(),
        help = 'Generate in silico restriction digest of reference FASTA.', 
        parents = [base_parser, txt_parser],
        formatter_class = formatter_class,
        epilog = epilog)
    digest_parser.add_argument(
        'infile', nargs = '?', default = '-',
        help = 'Input reference in FASTA format.')
    requiredNamed_digest = digest_parser.add_argument_group(
        'required named arguments')
    requiredNamed_digest.add_argument(
        '-r', '--restriction', required = True,
        type = restriction_seq, 
        help = '''Restriction cut sequence with "^" to indicate cut site.
                  e.g. Mbol = ^GATC''')
    digest_parser.set_defaults(function = hictools_digest.digest)
    commands[digest_command] = digest_parser
    
    # Truncate sub-parser
    truncate_command = 'truncate'
    truncate_parser = subparsers.add_parser('truncate',
        description = hictools_truncate.description(),
        help = 'Truncate FASTQ sequences at restriction enzyme ligation site.', 
        parents = [base_parser, txt_parser],
        formatter_class = formatter_class,
        epilog = epilog)
    truncate_parser.add_argument(
        'infile', nargs = '?', default = '-',
        help = 'Input file in FASTQ format.')
    truncate_parser.add_argument(
        '-s', '--summary', nargs = '?', 
        default = f'truncation_summary_{time.strftime("%Y%m%d-%H%M%S")}.txt', 
        help = 'Truncation summary file name.')
    truncate_parser.add_argument(
        '-n', '--sample', default = None,
        help = 'Sample name for truncation summary file.')
    requiredNamed_truncate = truncate_parser.add_argument_group(
        'required named arguments')
    requiredNamed_truncate.add_argument(
        '-r', '--restriction', required = True,
        type = restriction_seq, 
        help = '''Restriction cut sequence with "^" to indicate cut site.
                  e.g. Mbol = ^GATC''')
    truncate_parser.set_defaults(function = hictools_truncate.truncate)
    commands[truncate_command] = truncate_parser
    
    # Process sub-parser
    process_command = 'process'
    process_parser = subparsers.add_parser(process_command,
        description = hictools_process.description(),
        help = 'Determine HiC fragment mappings from '
               'named-sorted SAM/BAM file.', 
        parents = [base_parser, sam_parser],
        formatter_class = formatter_class,
        epilog = epilog)
    process_parser.add_argument(
        '-S', '--sam', dest = 'sam_out',
        action = 'store_true',
        help = 'Output alignments in SAM format.')
    process_parser.add_argument(
        '-u', '--gunzip', 
        action = 'store_true', dest = 'read_gzip',
        help = 'Read gzip compressed digest file.')
    requiredNamed_process = process_parser.add_argument_group(
        'required named arguments')
    requiredNamed_process.add_argument(
        '-d', '--digest', required = True,
        help = 'Output of hictools digest using same '
               'reference genome as used to map reads.')
    process_parser.set_defaults(function = hictools_process.process)
    commands[process_command] = process_parser
                  
    # Extract sub-parser
    extract_command = 'extract'
    extract_parser = subparsers.add_parser(extract_command,
        description = hictools_extract.description(),
        help = 'Extract HiC information encoded by hic process from SAM/BAM.', 
        parents = [base_parser, sam_parser],
        formatter_class = formatter_class,
        epilog = epilog)
    extract_parser.add_argument(
        '-z', '--gzip', 
        action = 'store_true', dest = 'write_gzip',
        help = 'Compress output using gzip')
    extract_parser.add_argument(
        '-n', '--sample', default = None,
        help = 'Sample name for input.')
    extract_parser.set_defaults(function = hictools_extract.extract)
    commands[extract_command] = extract_parser
    
    # Filter sub-parser
    filter_command = 'filter'
    filter_parser = subparsers.add_parser(filter_command,
        description = hictools_filter.description(),
        help = 'Filter SAM/BAM file processed with hictools process.', 
        parents = [base_parser, sam_parser],
        formatter_class = formatter_class,
        epilog = epilog)
    filter_parser.add_argument(
        '--min_inward', default = None, 
        type = positive_int, 
        help = 'Specify mininum insert size for inward facing read pairs.')
    filter_parser.add_argument(
        '--min_outward', default = None, 
        type = positive_int, 
        help = 'Specify mininum insert size for outward facing read pairs.')
    filter_parser.add_argument(
        '--max_ditag', default = None, 
        type = positive_int, 
        help = 'Specify maximum ditag size for read pairs.')
    filter_parser.add_argument(
        '-S', '--sam', dest = 'sam_out',
        action = 'store_true',
        help = 'Output alignments in SAM format.')
    filter_parser.set_defaults(function = hictools_filter.filter)
    commands[filter_command] = filter_parser
                  
    args = parser.parse_args()

    try:
        func = args.function
    except AttributeError:
        parser.print_help()
        sys.exit()
    
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
        # Print command specific sub-parser help.
        commands[args.command].print_help()
        sys.exit(1)

    args_dict = vars(args)
    [args_dict.pop(key) for key in ['command', 'function', 'verbose', 'log']]
    func(**vars(args))

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

def positive_int(value):
    ''' Custom argparse type for positive integer. '''
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(
            f'{value} is not a positive integer.')
    return ivalue

if __name__ == '__main__':
    try:
        main()
    finally:
        logging.shutdown()

