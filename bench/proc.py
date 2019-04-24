#!/usr/bin/env python3
import glob, sys, numpy as np

data = {}

for fname in sys.argv[1:]:
    with open(fname, "r") as f:
        lines = f.readlines()
        rows = [[int(x) for x in line.split(" ")] for line in lines]
        data[fname] = np.array(rows)

#import pdb; pdb.set_trace()

# rank the keys by best time
sumtimes = [(k,sum(rows.T[-1])) for k, rows in data.items()]

sumtimes.sort(key=lambda x: x[1])

print(sumtimes)
