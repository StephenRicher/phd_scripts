#!/usr/bin/env python3


''' Generate an in-silico restriction enzyme digest of a genome in 
FASTA format. 
'''


import argparse, sys, logging, re, select
from functools import partial

sys.path.insert(0, '/home/stephen/phd/scripts/python3_modules/')
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
        help = '''Restriction cut sequence with "^" to indicate cut site.
                  e.g. Mbol = ^GATC''')

    args = parser.parse_args()

    data_in_stdin = select.select([sys.stdin,],[],[],0.0)[0]
    if not data_in_stdin and args.infile == '-':
        sys.stderr.write(f'Error: No input provided.\n\n')
        parser.print_help()
        sys.exit(1)

    return vars(args)

def find_cut_sites(ref_seq, ref, restriction, out_obj):
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')

    if not ref_seq:
        log.error(f'Reference {ref} contains no sequence.')
        ec = 1
    else:
        overhang = restriction.index('^')
        site = restriction.replace('^','')
        matches = re.finditer(site, ref_seq)
        index = 0
        previous_end = 0
        for match in matches:
            # Skip if restriction sequence at start of reference.
            if match.start() == 0:
                continue
            index += 1
            start = 1 if index == 1 else previous_end + 1
            end = match.start() + overhang
            out_obj.write(f'{ref}\t{start}\t{end}\t{index}\n')    
            previous_end = end
        out_obj.write(
            f'{ref}\t{previous_end + 1}\t{len(ref_seq)}\t{index + 1}\n')
        ec = 0
    return ec
    
def check_valid_seq(seq):
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')
    
    return re.search('[^ATCGURYKMSWBDHVN-]', seq.strip('\n'), re.IGNORECASE)
        
def main(infile, output, restriction):
    
    ''' Iterate through each infile. '''
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')

    with smart_open(output, 'w') as out_obj, \
            smart_open(in_path, 'r') as in_obj:
        log.info(f'Writing output to {output}.')
        header = 1
        for index, line in enumerate(in_obj):
            if line.startswith('>'):
                if header > 1:
                    find_cut_sites(''.join(seqs), ref, restriction, out_obj)
                ref = line.rsplit()[0][1:]
                log.info(f'Digesting reference {ref}.')
                header += 1
                seqs = []
            elif header == 1 and index == 0:
                log.error(f'FASTA line 1 does not begin with ">".')
                sys.exit(1)
            else:
                invalid = check_valid_seq(line)
                if invalid:
                    log.error(f'Invalid FASTA character {invalid.group(0)} '
                               'on line {index}.')
                    sys.exit(1)
                seqs.append(line.upper().strip('\n'))
        find_cut_sites(''.join(seqs), ref, restriction, out_obj)
                    
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
        
    finally:
        logging.shutdown()

