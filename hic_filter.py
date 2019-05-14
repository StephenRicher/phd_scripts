#!/usr/bin/env python3

import collections
import bisect
import time
import re

class sam:
    
    def __init__(self, sam_line):
        self.qname = sam_line[0]
        self.flag = int(sam_line[1])
        self.rname = sam_line[2]
        self.left_pos = int(sam_line[3])
        self.mapq = int(sam_line[4])
        self.cigar = sam_line[5]
        self.rnext = sam_line[6]
        self.pnext = int(sam_line[7])
        self.tlen = int(sam_line[8])
        self.seq = sam_line[9]
        self.qual = sam_line[10]
        self.optional = sam_line[11:]
        self.fragment = None
    
    @property    
    def is_reverse(self):
        return True if (self.flag & 0x10 != 0) else False
    
    @property
    def is_read1(self):
        return True if (self.flag & 0x40 != 0) else False
    
    @property
    def is_paired(self):
        return True if (self.flag & 0x1 != 0) else False
    
    @property
    def reference_length(self):
        cigar_split = re.split("(\d+)", self.cigar)[1:]
        length = 0
        for idx, val in enumerate(cigar_split):
            if idx & 1 and val not in ["I", "S", "H", "P"]:
                length += int(cigar_split[idx-1])
        return length
    
    @property
    def right_pos(self):
        return self.left_pos + (self.reference_length - 1)
    
    @property
    def five_prime_pos(self):
        if self.is_reverse:
            return self.right_pos
        else:
            return self.left_pos
    
    @property
    def three_prime_pos(self):
        if self.is_reverse:
            return self.left_pos
        else:
            return self.right_pos
    
    @property
    def middle_pos(self):
        return round((self.right_pos + self.left_pos)/2)
    
    def print_sam(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            self.qname, self.flag, self.rname, self.left_pos, self.mapq, self.cigar,
            self.rnext, self.pnext, self.tlen, self.seq, self.qual, '\t'.join(map(str,self.optional)))

class fragment:
    def __init__(self, fragment_number, start, end):
        self.fragment_number = fragment_number
        self.start = start
        self.end = end

def is_valid(read1, read2):
    if read1.qname != read2.qname:
        print ("Error: Qname mismatch: {} {} - is file name sorted?".format(read1.qname, read2.qname))
    elif not read1.is_paired and not read2.is_paired:
        print("Error: {} is not paired".format(read1.qname))
    elif read1.is_read1 == read2.is_read1:
        print("Error: R1 and R2 flags in {} not correctly set".format(read1.qname))
    elif read1.pnext != read2.left_pos or read2.pnext != read1.left_pos:
        print("Error: Mate position mismatch in {}.".format(read1.qname))
    else:
        return True
    return False
  
def is_trans(read1, read2):
    """ Return true if reads do not align to same reference """
    return read1.rname != read2.rname
    
def get_fragment(read, digest):
    rf_num = bisect.bisect_left(digest[read.rname], read.middle_pos)
    rf_start = 1 if rf_num == 0 else digest[read.rname][rf_num - 1] + 1
    rf_end = digest[read.rname][rf_num]
    read.fragment = rf_num
    return fragment(rf_num, rf_start, rf_end)
    
def tag_length(read, fragment):
    if read.is_reverse:
        return read.five_prime_pos - fragment.start + 1
    else:
        return fragment.end - read.five_prime_pos + 1

def process_digest(digest_filename):
    with open(digest_filename) as digest:
        d = {}
        for fragment in digest:
            [chromosome, start, end, fragment_number] = fragment.split()
            if chromosome not in d.keys():
                d[chromosome] = []
            d[chromosome].append(int(end))
    return(d)

def reorder_read_pair(read1, read2):
    """
    Return a pair of reads such that read1 is left of read 2.
    Read pairs aligning to different chromosomes are returned unchanged.
    """
    
    if not is_trans(read1, read2) and read1.left_pos > read2.left_pos:
        r1_reorder = read2
        r2_reorder = read1
    else:
        r1_reorder = read1
        r2_reorder = read2
    return r1_reorder, r2_reorder

def get_orientation(read1, read2):
    """
    Return relative orientation of read pairs.
    Assumes read pairs have been ordered such that read 1 is five prime of read 2.
    """
    
    if read1.is_reverse:
        if read2.is_reverse:
            orientation = "same_reverse"
        else:
            orientation = "outward"
    else:
        if read2.is_reverse:
            orientation = "inward"
        else:
            orientation = "same_forward"
    return orientation
       
def on_same_fragment(read1, read2):
    if read1.fragment == None or read2.fragment == None:
        sys.exit("Error read fragment not set for {}.".format(read1.qname))
    return read1.fragment == read2.fragment

def on_adjacent_fragments(read1, read2):
    if read1.fragment == None or read2.fragment == None:
        sys.exit("Error read fragment not set for {}.".format(read1.qname))
    return abs(read1.fragment - read2.fragment) == 1

def run_filter(read1, read2, digest):
    if not is_valid(read1, read2):
        status = "invalid_format"
    else:
        read1, read2 = reorder_read_pair(read1, read2)
        orientation = get_orientation(read1, read2)
        read1_fragment = get_fragment(read1, digest)
        read2_fragment = get_fragment(read2, digest)
        ditag_length = tag_length(read1, read1_fragment) + tag_length(read2, read2_fragment)
        if ditag_length < 150 or ditag_length > 800:
            status = "wrong_size"
        elif is_trans(read1, read2):
            status = "valid_trans"
        elif on_same_fragment(read1, read2):
            if orientation == "outward":
                status = "self_circle"
            elif orientation == "inward":
                status = "dangling_end"
            else:
                status = "same_strand_error"
        elif orientation == "inward" and on_adjacent_fragments(read1, read2):
            status = "religation"
        elif orientation == "inward" and 150 < read2.five_prime_pos - read1.five_prime_pos + 1 < 800:
            status = "contiguous"
        else:
            status = "valid_cis"
        return status
                

def main():
    digest_filename="/home/stephen/x_am/RC-BB1219/stephen_test/genomes/GRCh38/digest/GRCh38_primary_assembly_Mbo1_digest.txt"
    d = process_digest(digest_filename)
    stats = {"total_pairs": 0, "invalid_format": 0, "wrong_size": 0, "self_circle": 0, "dangling_end": 0, 
        "same_strand_error": 0, "religation": 0, "contiguous": 0, "valid_cis": 0, "valid_trans": 0}
        
    with open("/media/stephen/Data/hic_analysis/data/HB2_CL4-1.fix.pair.sam") as sam_file:
        for line in sam_file:
            if line.startswith("@"):
                continue
            else:
                try:
                    read1 = sam(line.split())
                    read2 = sam(next(sam_file).split())
                except StopIteration:
                    print ("Error: odd number of fastq sequences in file")
                status = run_filter(read1, read2, d)
                
                read1.optional.append("hf:Z:{}".format(status))
                read2.optional.append("hf:Z:{}".format(status))
                stats["total_pairs"] += 1
                stats[status] += 1
                print(read1.print_sam())
                print(read1.print_sam())
        for key,val in stats.items():
            print (key, "\t", val)

def pysam_test():
    import pysam
    samfile = pysam.AlignmentFile("/media/stephen/Data/hic_analysis/data/test.sam", "rb")
    mysamfile = open("/media/stephen/Data/hic_analysis/data/test2.sam")
    for read_pysam, read_mysam in zip(samfile.fetch(), mysamfile):
        myread = sam(read_mysam.split())
        assert read_pysam.get_reference_positions()[-1] + 1 == myread.right_pos
        assert read_pysam.get_reference_positions()[0] + 1 == myread.left_pos
        assert read_pysam.reference_length == myread.reference_length
        assert read_pysam.is_paired == myread.is_paired
        assert read_pysam.is_reverse == myread.is_reverse
        if myread.is_reverse:
            assert read_pysam.get_reference_positions()[-1] + 1 == myread.five_prime_pos
            assert read_pysam.get_reference_positions()[0] + 1 == myread.three_prime_pos
        else:
            assert read_pysam.get_reference_positions()[-1] + 1 == myread.three_prime_pos
            assert read_pysam.get_reference_positions()[0] + 1 == myread.five_prime_pos
    mysamfile.close()
            
if __name__ == "__main__":
    main()
    #pysam_test()
    

