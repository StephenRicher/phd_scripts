#!/usr/bin/env python3

import re, argparse, sys, gzip, contextlib, os, stat, select, logging

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', nargs = '?', help = 'Specify input file name.', type = argparse.FileType('r'), default = sys.stdin)
parser.add_argument('-l', '--log', nargs = '?', help = 'Specify log file name', type = argparse.FileType('w'),  default = sys.stderr)
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

def run():
    ''' Run main loop. '''
    try:
        if not args.file and sys.stdin.isatty():
            parser.print_help()

        invalid = False

        with args.file as f:
            for index, line in enumerate(f):
                line = line.rstrip('\n')
                in_header = True if line.startswith('>') else False
                header_above = None
                if not line:
                    log.warning(f'Error: Line {index} is a blank line.')
                    invalid = True
                elif not in_header:
                    invalid_characters = re.search('[^ATCGURYKMSWBDHVN-]', line, re.IGNORECASE)
                    if invalid_characters:
                        log.warning(f'Error: Invalid sequence character "{invalid_characters.group(0)}" at line {index}.')
                        invalid = True
                    elif in_header == header_above:
                        log.warning('Error: FASTA sequence contains line wrapping.')
                        invalid = True
                # If we are in header and the previous line was also a header...
                elif in_header == header_above:
                    log.warning(f'Error: Header missing sequence at line {index}.')
                    invalid = True
                header_above = in_header
            if in_header:
                log.warning('Error: Header missing sequence on last last.')
                invalid = True
    finally:
        logging.shutdown()
        if invalid:
            sys.exit(1)
                
if __name__ == '__main__':
    run()
