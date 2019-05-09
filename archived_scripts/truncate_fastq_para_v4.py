import gzip
import argparse
import sys
import os
import multiprocessing
import time
import contextlib
import datetime

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", help = "Input fastq file(s)", nargs = "*")
parser.add_argument("-r", "--restriction_sequence",     help = "Specify restriction sequence and cut site. e.g. '^GATC' for Mbol", type = str)
parser.add_argument('--gz', help = "Specify if stdin is gzipped.", action='store_true')
parser.add_argument('--stdout', help = "Specify if output should be sent to stdout.", action='store_true')
args = parser.parse_args()

if not args.infile:
        sys.exit("Please provide an input file.")
if not args.restriction_sequence:
        sys.exit("Please provide restriction sequence.")

# NEED TO PUT THIS INTO A FUCTION PARAMETER!
summary_file = "truncate_summary_{0:%Y-%m-%d_%H-%M-%S}.txt".format(datetime.datetime.now())

@contextlib.contextmanager
def smart_open (filename, gz = False, mode = None):
    """ Custom context file opener which allows which can read text and gz compressed data
        from files and stdin. Also can write text to stdout write gz compressed data to file."""
    if mode not in ["w", "r"]:
        sys.exit("Smart open mode not set - define read ('r') or write ('w').")
    if filename and filename != '-':
        if mode == "w" or is_gz_file(filename):
            fh = gzip.open(filename, mode + 't')
        else:
            fh = open(path, mode + 't')
    else:
        if mode == "r":
            if not sys.stdin.isatty():
                if gz:
                    fh = gzip.open(sys.stdin.buffer, mode + 't')
                else:
                    fh = sys.stdin
            else:
                sys.exit("Nothing in stdin")
        else:
            fh = sys.stdout
            
    try:
        yield fh
    finally:
        if fh not in [sys.stdout, sys.stdin]:
            fh.close()

def is_gz_file (filepath):
    """ Function doc """
    try:
        with open(filepath, 'rb') as f:
            return f.read(2).hex() == '1f8b' and filepath.endswith("gz")
    except IOError as e:
            sys.exit("Error: " + str(e))

def process_fastq(lines = None):
    """ Reads 4 line list of fastq record into dictionary. """
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}

def truncate_fastq (record, ligation_sequence, truncation_index):
    """ Function doc """
    is_truncated = False
    if ligation_sequence in record['sequence'].upper():
        truncation_site = record['sequence'].upper().index(ligation_sequence) + truncation_index
        record['sequence'] = record['sequence'][0:truncation_site]
        record['quality'] = record['quality'][0:truncation_site]
        is_truncated = True
    return (record, is_truncated)

def process_restriction (restriction_sequence):
    """ Function doc """
    cut_position = restriction_sequence.index("^")
    sequence = restriction_sequence.replace("^", "").upper()
    ligation_sequence = sequence + sequence[cut_position:]
    truncation_index = len(sequence)
    return ligation_sequence, truncation_index
    
def iterate_fastq (f_in, outfile, restriction_sequence):
    """ Function doc """
    try:
        with smart_open(filename = outfile, gz = True, mode = "w") as f_out:
            ligation_sequence, truncation_index = process_restriction(restriction_sequence)
            total_truncated = 0
            total_len = 0
            total_reads = 0
            lines = []
            try:
                for line in f_in:
                    lines.append(line.rstrip())
                    if len(lines) == 4:
                        record, is_truncated = truncate_fastq(process_fastq(lines), ligation_sequence, truncation_index)
                        lines = []
                        total_reads += 1
                        if is_truncated:
                            total_truncated += 1
                            total_len = total_len + len(record['sequence']) - 1
                        f_out.write('{0}\n{1}\n{2}\n{3}\n'.format(record['name'], record['sequence'], record['optional'], record['quality']))
            except Exception as e:
                sys.exit("Error: " + str(e) + " - " + f_in.name)
        try:
            with open(summary_file, "a+") as summary:
                if os.stat(summary_file).st_size == 0:
                    summary.write("file" + '\t' + "total_reads_processed" + '\t' + "truncated" + '\t' + "mean_truncated_length" + "\n")
                summary.write(f_in.name + '\t' + str(total_reads) + '\t' + str(total_truncated) + '\t' + str(total_reads / total_len) + "\n")
        except IOError as e:
            sys.exit("Error: " + str(e) + " - " + summary_file)
    except IOError as e:
        sys.exit("Error: " + str(e) + " - " +  outfile)

def run_truncater (infile):
    """ Function doc """
    try:
        with smart_open(filename = infile, gz = args.gz, mode = "r") as f_in:
            if args.stdout:
                outfile = "-"
            elif infile == "-":
                outfile = "fastq_truncated.gz"
            else:
                outfile = os.path.splitext(f_in.name)[0] + '.trunc.gz'
            iterate_fastq(f_in, outfile, args.restriction_sequence)
    except IOError as e:
        sys.exit("Error: " + str(e))

if __name__ == '__main__':
    pool = multiprocessing.Pool(multiprocessing.cpu_count()-1 or 1)
    results = [pool.apply_async(run_truncater, args=(infile,)) for infile in args.infile]
    output = [p.get() for p in results]
    print(output)
