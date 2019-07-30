#!/usr/bin/env python3

import re

from class_validators import *

class sam:

    qname = RegexMatch(re.compile(
        r"^[!-?A-~]+$"))
    flag = IntRange(0, (2**16) - 1)
    rname = RegexMatch(re.compile(
        r"^\*$|^[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*$"))
    left_pos = IntRange(0, (2**31) - 1)
    mapq = IntRange(0, (2**8) - 1)
    cigar = RegexMatch(re.compile(
        r"^\*$|^([0-9]+[MIDNSHPX=])+$"))
    rnext = RegexMatch(re.compile(
        r"^\*$|^=$|^[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*$"))
    pnext = IntRange(0, (2**31) - 1)
    tlen = IntRange((-2**31) + 1, (2**31) - 1)
    seq = RegexMatch(re.compile(
        r"^\*$|^[A-Za-z=.]+$"))
    qual = RegexMatch(re.compile(
        r"^[!-~]+$"))
    
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
        self.optional = self.read_opt(sam_line[11:])
        
    def read_opt(self, all_opts):
        """ Process optional SAM files into a dictionary """
        d={}
        for opt in all_opts:
            tag_and_type = opt[0:opt.rindex(':')]
            type_ = opt.split(':')[1]
            value = opt[opt.rindex(':') + 1:]
            if type_ == 'i':
                value = int(value)
            elif type_ == 'f':
                value = float(value)
            d[tag_and_type] = value
        return d
        
    def print_opt(self, opt):
        """ Output optional SAM fields as tab-delimated string """
        opt_out = ""
        for tag_and_type, value in opt.items():
            opt_out += f'{tag_and_type}:{value}\t'
        return opt_out
    
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
        return self.left_pos if self.is_reverse else self.right_pos
    
    @property
    def middle_pos(self):
        return round((self.right_pos + self.left_pos)/2)
    
    def print_sam(self):
        return f'{self.qname}\t{self.flag}\t{self.rname}\t{self.left_pos}\t{self.mapq}\t{self.cigar}\t{self.rnext}\t{self.pnext}\t{self.tlen}\t{self.seq}\t{self.qual}\t{self.print_opt(self.optional)}\n'
