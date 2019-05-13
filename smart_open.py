#!/usr/bin/env python

import argparse, sys, gzip, contextlib, os, stat, select

@contextlib.contextmanager
def smart_open (filename, mode = "read"):
    """ Custom context file opener which can read from uncompressed and gzip compressed data from stdin (default)
        input files. Named pipes are auto detected and read but cannot be in gzip format. It can write 
        to stdout (default) or to a specified outfile. Outfiles suffixed with '.gz' are automatically gzipped. """
    
    def named_pipe (filename):
        """ Return true if input file is a named pipe """
        try:
            return stat.S_ISFIFO(os.stat(filename).st_mode)
        except OSError as e:
            sys.exit(str(e))
            
    def is_gzip (filename):
        """ Auto detects whether an input file is gz compressed.
            Checks for magic number and if filepaths ends in '.gz'"""
        try:
            f = open(filename, 'rb')
        except IOError as e:
            sys.exit(str(e))
        else:
            with f:
                return f.read(2).hex() == '1f8b' and filename.endswith("gz")
                
    if mode == "read":
        if filename and filename != '-':
            # If filename is a named pipe (i.e. via process substitution) 
            # then do not check if gzip because this will close the pipe. 
            # Otherwise check if file is gzipped.
            if not named_pipe(filename) and is_gzip(filename):
                try:
                    fh = gzip.open(filename, 'rt')
                except IOError as e:
                    sys.exit(str(e))
            else:
                try:
                    fh = open(filename, 'rt')
                except IOError as e:
                    sys.exit(str(e))
        else:
            if not sys.stdin.isatty():
                if sys.stdin.buffer.peek(2).hex()[0:4] == '1f8b':
                    fh = gzip.open(sys.stdin.buffer, mode = 'rt')
                else:
                    fh = sys.stdin
            else:
                sys.exit("Nothing in stdin and no infile provided.")
    elif mode == "write":
        if filename and filename != '-':
            if filename.endswith("gz"):
                try:
                    fh = gzip.open(filename, 'wt')
                except IOError as e:
                    sys.exit(str(e))
            else:
                try:
                    fh = open(filename, 'wt')
                except IOError as e:
                    sys.exit(str(e))
        else:
            fh = sys.stdout
    else:
        sys.exit("Specify mode as either 'read' or 'write'.")
            
    try:
        yield fh
    finally:
        if fh not in [sys.stdin, sys.stdout]:
            fh.close()

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help = "Specify input file name", nargs = "?")
    parser.add_argument("-o", "--out", help = "Specify output file name", nargs = "?")
    args = parser.parse_args()
    
    if not args.file and sys.stdin.isatty():
        parser.print_help()

    # Test functionality
    with smart_open(args.file, mode = "read") as f, smart_open(args.out, mode = "write") as out:
        try:
            for line in f:
                out.write(line)
        except UnicodeError as e:
            sys.exit("Error: {}\nPossibly compressed data in process substitution.".format(str(e)))
            
