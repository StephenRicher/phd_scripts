#!/usr/bin/env python3

import sys, pytest
sys.path.append("/home/stephen/h/phd/scripts2/")

from hic_filter import *

# Upload test digest file
@pytest.fixture
def digest_file():
    return process_digest("/home/stephen/h/phd/scripts2/testing/hic_filter/test_digest.txt")

# Read test dataset into two lists to be passed to parameters.
with open("/home/stephen/h/phd/scripts2/testing/hic_filter/test_selfcircle.sam") as f:
    reads1 = []
    reads2 = []
    for line in f:
        if line.startswith("@"):
            continue
        else:
            reads1.append(sam(line.split()))
            reads2.append(sam(next(f).split()))

# Iterate through each pair of parameters
@pytest.mark.parametrize('read1', reads1)
@pytest.mark.parametrize('read2', reads2)
class TestParametrized:

    @pytest.mark.parametrize('output', [True])
    def test_is_valid(self, read1, read2, output):
        assert is_valid(read1, read2) == output

    @pytest.mark.parametrize('output', ['self_circle'])
    def test_status(self, read1, read2, output, digest_file):
        assert run_filter(read1, read2, digest_file) == output
