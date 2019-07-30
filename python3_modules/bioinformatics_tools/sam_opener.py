#!/usr/bin/env python3

""" Custom context manager for reading and writing SAM/BAM files.
    Requires samtools to be installed.
"""

import sys, logging, contextlib, subprocess

@contextlib.contextmanager
def sam_open(filename: str = None, mode: str = 'r', samtools = 'samtools'):
    
    """ Custom context manager for reading and writing SAM/BAM files. 
    """
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')
    
    if mode not in ['r', 'w', 'wt', 'wb', 'rt', 'rb']:
        log.error(f'Invalid mode {mode} for sam_open.')
    
    if 'r' in mode:
        try:
            if filename == '-':
                stdin = sys.stdin
            else:
                stdin = None
            p = subprocess.Popen(
                ['samtools', 'view', '-h', filename], 
                stdout = subprocess.PIPE, 
                stdin = stdin, encoding = 'utf8')
            fh = p.stdout

        except FileNotFoundError:
            log.error(f'Samtools not found.')
            sys.exit(1)
    else:
        try:
            if filename == '-':
                outfile = sys.stdout.buffer
            else:
                try:
                    outfile = open(filename, 'wt')
                except IOError:
                    log.exception(f'Unable to open {filename}.')
                    sys.exit(1)
            out = '-b' if 'b' in mode else '-h'
            p = subprocess.Popen(
                    ['samtools', 'view', out], stdout = outfile,
                    stdin = subprocess.PIPE, encoding = 'utf8')
            if filename != '-':
                outfile.close()
            fh = p.stdin
        except FileNotFoundError:
            log.error(f'Samtools not found.')

    try:
        yield fh
    finally:
        try:
            fh.close()
        except AttributeError:
            pass
