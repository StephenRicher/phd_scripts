#!/usr/bin/env python3

import math, sys, bisect

from class_validators import *

class fragment:
    
    fragment_number = IntRange(0, math.inf)
    fragment_start = IntRange(0, math.inf)
    fragment_end = IntRange(0, math.inf)
    
    def __init__(self, fragment_number, start, end):
        self.fragment_number = fragment_number
        self.start = start
        self.end = end

def is_valid(read1, read2):
    if read1.qname != read2.qname:
        sys.stderr.write(f'Qname mismatch: {read1.qname} {read2.qname} - is file name sorted?')
    elif not read1.is_paired and not read2.is_paired:
        sys.stderr.write(f'{read1.qname} is not paired')
    elif read1.is_read1 == read2.is_read1:
        sys.stderr.write(f'R1 and R2 flags in {read1.qname} not correctly set')
    elif read1.pnext != read2.left_pos or read2.pnext != read1.left_pos:
        sys.stderr.write(f'Mate position mismatch in {read1.qname}.')
    else:
        return True
    return False
    
def interaction_type(read1, read2):
    if read1.rname != read2.rname:
        interaction = "trans"
    else:
        interaction = "cis"
    return interaction
    
def get_fragment(read, digest):
    rf_num = bisect.bisect_left(digest[read.rname], read.middle_pos)
    rf_start = 1 if rf_num == 0 else digest[read.rname][rf_num - 1] + 1
    rf_end = digest[read.rname][rf_num]
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
            [ref, start, end, fragment_number] = fragment.split()
            if not (int(start) > 0 and int(end) > 0):
                raise ValueError(f'Negative fragment start/end positions on ref {ref}.')
            if ref not in d.keys():
                if not (int(start) == 1 and int(fragment_number) == 1):
                    raise ValueError(f'Invalid first fragment in ref {ref}.')
                d[ref] = []
            d[ref].append(int(end))
    return(d)

def reorder_read_pair(read1, read2):
    """
    Return a pair of reads such that read1 is left of read 2.
    Read pairs aligning to different chromosomes are returned unchanged.
    """
    
    if interaction_type(read1, read2) == "cis" and read1.left_pos > read2.left_pos:
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
            orientation = "Same-reverse"
        else:
            orientation = "Outward"
    else:
        if read2.is_reverse:
            orientation = "Inward"
        else:
            orientation = "Same-forward"
    return orientation
       
def on_same_fragment(r1_fragment, r2_fragment):
    return r1_fragment.fragment_number == r2_fragment.fragment_number

def on_adjacent_fragments(r1_fragment, r2_fragment):
    return abs(r1_fragment.fragment_number == r2_fragment.fragment_number) == 1

def run_filter(read1, read2, digest):
    if not is_valid(read1, read2):
       raise ValueError(f'Invalid format in {read1.qname}.')     
    read1, read2 = reorder_read_pair(read1, read2)
    orientation = get_orientation(read1, read2)
    read1_fragment = get_fragment(read1, digest)
    read2_fragment = get_fragment(read2, digest)
    ditag_length = tag_length(read1, read1_fragment) + tag_length(read2, read2_fragment)
    insert_size = read2.right_pos - read1.left_pos + 1
    interaction = interaction_type(read1, read2)
    fragment_seperation = abs(read2_fragment.fragment_number - read1_fragment.fragment_number)
    # ADD this to a dictionary??
    return orientation, ditag_length, insert_size, interaction, fragment_seperation
