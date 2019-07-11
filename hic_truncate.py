#!/usr/bin/env python3

import argparse, sys, re

def restriction_seq(value):
    ''' Custom argparse type restriction enzyme type. '''
    if value.count('^') != 1:
        raise argparse.ArgumentTypeError(f'Error: Restriction site {value} must contain one "^" to indicate cut site.')
    elif re.search('[^ATCG^]', value, re.IGNORECASE)
        raise argparse.ArgumentTypeError(f'Error: Restriction site {value} may only contain "ATCG^".')
    else:
        return value

parser = argparse.ArgumentParser(description = 'Script to truncate proximity ligated restriction fragments in FASTQ files.')
parser.add_argument('-f', '--file', nargs = '?', help = 'Input file', type = argparse.FileType('r'), default = sys.stdin)
parser.add_argument('-o', '--out', nargs = '?', help = 'Output file', type = argparse.FileType('w'),  default = sys.stdout)
parser.add_argument('-l', '--log', nargs = '?', help = 'Log file', type = argparse.FileType('w'),  default = sys.stderr)
parser.add_argument('-r', '--results', nargs = '?', help = 'Results file', type = argparse.FileType('w'),  default = sys.stderr)
parser.add_argument('-s', '--sample', help = 'Specify sample name.', type = str, default = None)
# Create group for required named arguments.
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-e', '--restrictionSite', type = restriction_seq, help = 'Restriction cut site (e.g. Mbo1: ^GATC).', required = True)
args = parser.parse_args()

def handle_exception(exc_type, exc_value, exc_traceback):
    ''' Handle uncaught exceptions (excluding KeyboardInterrupt) 
        through the logging module. '''
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logging.getLogger(__name__).critical("Uncaught exception", exc_info = (exc_type, exc_value, exc_traceback))

sys.excepthook = handle_exception

logging.basicConfig(stream = args.log, format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s', level = logging.DEBUG)
log = logging.getLogger(__name__)

def process_restriction(restriction):
    restriction_seq = restriction.replace('^','')
    cut_site1 = restriction.index('^')
    ligation2 = restriction_seq[cut_site1:]
    cut_site2 = len(restriction) - cut_site1 - 1
    ligation1 = restriction_seq[0:cut_site2]
    ligation_seq = ligation1 + ligation2
    cut_index = len(ligation1)
    return ligation_seq, restriction_seq

def run():
    ''' Run main loop. '''
    try:
        if not args.file and sys.stdin.isatty():
            parser.print_help()
 
        if not args.sample:
            if args.file is not sys.stdin:
                args.sample = args.file
            else:
                args.sample = 'sample'
 
        ligation_seq, restriction_seq = process_restriction(args.restrictionSite)
        truncated = 0
        truncated_length = 0
    
        with args.file as f, args.out as out, args.results as results:
            for index, line in enumerate(f):
                line=line.rstrip('\n')
                # Sequence line
                if index % 4 == 1 and ligation_seq in line:
                    line = line[0: line.index(ligation_seq)] + restriction_seq
                    is_truncated = True
                seq_length = len(line)
                # Quality line
                elif index % 4 == 3:
                    line = line[0:seq_length]
                    if truncated:
                        truncated += 1
                        truncated_length += seq_length
                        is_truncated = False
                out.write(f'{line}\n')
            results.write(f'{args.sample}\t{truncated}\t{truncated_length/truncated}\n')
    finally:
        logging.shutdown()

if __name__ == '__main__':
    run()
    
            
