#!/usr/bin/env python3

import sys, contextlib

@contextlib.contextmanager
def smart_read (filename):
    """ Custom context file opener which can read from from stdin or input file."""
    if filename != sys.stdin:
        try:
            fh = open(filename, 'rt')
        except IOError as e:
            sys.exit(str(e))
    else:
        fh = sys.stdin
    try:
        yield fh
    finally:
        if fh is not sys.stdin:
            fh.close()
