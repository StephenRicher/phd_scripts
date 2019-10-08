#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd

def main():
    
    assert len(sys.argv) == 3

    sample1 = sys.argv[1]
    sample2 = sys.argv[2]
    
    with open(sample1) as f:
        track = f.readline()
        
    a = pd.read_csv(sample1, delimiter = '\t', skiprows = 1, 
        header = None).set_index([0,1,2])
    b = pd.read_csv(sample2, delimiter = '\t', skiprows = 1, 
        header = None).set_index([0,1,2])
                
    diff = a.sub(b).dropna()
    
    output = track + diff.to_csv(sep = '\t', header = False)
    
    sys.stdout.write(output)

if __name__ == "__main__":
    main()
