#!/usr/bin/env python3

"""
Custom context managers for reading and writing GZIP files.
File objects can be reopened for GZIP read/writing
"""

import io, os, sys, logging, gzip, contextlib, binascii, stat, subprocess

def is_gzip(filepath):
    
    ''' Check for GZIP magic number byte header. '''
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')
    
    with open(filepath, 'rb') as f:
        return binascii.hexlify(f.read(2)) == b'1f8b'
        
def named_pipe(path):
    
    """ Check if file is a named pipe. """
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')
    
    if stat.S_ISFIFO(os.stat(path).st_mode):
        pipe = True
    else:
        pipe = False
    return pipe

@contextlib.contextmanager
def smart_open(
    filename: str = None, mode: str = 'r', gz = False, *args, **kwargs):
    
    """ Custom context manager for reading and writing uncompressed and
        GZ compressed files. 
        
        Interprets '-' as stdout or stdin as appropriate. Also can auto 
        detect GZ compressed files except for those inputted to stdin or 
        process substitution. Uses gzip via a subprocess to read/write
        significantly faster than the python implementation of gzip.
        On systems without gzip the method will default back to the python
        gzip library.
        
        Ref: https://stackoverflow.com/a/45735618
    """

    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')

    # If decompress not set then attempt to auto detect GZIP compression.
    if (    not gz and
            'r' in mode and
            filename.endswith('.gz') and 
            not named_pipe(filename) and 
            is_gzip(filename)):
        log.info(
            f'{filename} detected as gzipped. Decompressing...') 
        gz = True

    if gz:
        encoding = None if 'b' in mode else 'utf8'
        if 'r' in mode:
            try:
                if filename == '-':
                    stdin2 = sys.stdin
                else:
                    stdin2 = None
                p = subprocess.Popen(
                    ['zcat', filename], stdout = subprocess.PIPE, 
                    stdin = stdin2, encoding = encoding)
                fh = p.stdout
            except FileNotFoundError:
                if filename == '-':
                    infile = sys.stdin.buffer
                else:
                    infile = filename
                fh = gzip.open(infile, mode, *args, **kwargs)
        else:
            try:
                if filename == '-':
                    outfile = sys.stdout.buffer
                else:
                    try:
                        outfile = open(filename, mode)
                    except IOError:
                        log.exception(f'Unable to open {filename}.')
                        sys.exit(1)
                p = subprocess.Popen(
                        ['gzip', '-f'], stdout = outfile,
                        stdin = subprocess.PIPE, encoding = encoding)
                if filename != '-':
                    outfile.close()
                fh = p.stdin
            except FileNotFoundError:
                if filename == '-':
                    outfile = sys.stdout.buffer
                else:
                    outfile = filename
                fh = gzip.open(outfile, mode, *args, **kwargs)
    else:
        if filename == '-':
            if 'r' in mode:
                stream = sys.stdin
            else:
                stream = sys.stdout
            if 'b' in mode:
                fh = stream.buffer
            else:
                fh = stream
        else:
            fh = open(filename, mode, *args, **kwargs)

    try:
        yield fh
    finally:
        try:
            fh.close()
        except AttributeError:
            pass
