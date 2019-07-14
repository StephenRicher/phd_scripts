#!/usr/bin/env python3

"""
Custom context managers for reading and writing GZIP files.
File objects can be reopened for GZIP read/writing
"""

import os, sys, logging, gzip, contextlib, binascii, stat

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
def gzip_read(path = None, gunzip = False):
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')

    # If decompress not set then attempt to auto detect GZIP compression.
    if (    not gunzip and 
            path.endswith('.gz') and 
            not named_pipe(path) and 
            is_gzip(path)):
        log.info(f'Input {path} detected as gzipped. Decompressing...') 
        gunzip = True
        
    try:
        if gunzip:
            if path == '-':
                fobj = gzip.open(sys.stdin.buffer, mode = "rt")
            else:
                fobj = gzip.open(path, mode = "rt")
        elif path == '-':
            fobj = sys.stdin
        else:
            fobj = open(path, mode = "rt")
    except IOError:
        log.exception(f'Unable to open {path}.')
        sys.exit()
    try:
        yield fobj
    finally:
        if path is not '-':
            fobj.close()
            
@contextlib.contextmanager
def gzip_write(path = None, compress = False):
    
    fun_name = sys._getframe().f_code.co_name
    log = logging.getLogger(f'{__name__}.{fun_name}')

    try:
        if compress:
            if path == '-':
                fobj = gzip.open(sys.stdout.buffer, mode = "wt")
            else:
                fobj = gzip.open(path, mode = "wt")
        elif path == '-':
            fobj = sys.stdout
        else:
            fobj = open(path, mode = "rt")
    except IOError:
        log.exception(f'Unable to open {path}.')
        sys.exit()
    try:
        yield fobj
    finally:
        if path is not '-':
            fobj.close()
