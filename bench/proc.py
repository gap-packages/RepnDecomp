#!/usr/bin/env python3
import glob, sys, numpy as np

data = {}

for fname in sys.argv[1:]:
    with open(fname, "r") as f:
        lines = f.readlines()
        rows = [[int(x) for x in line.split(" ")] for line in lines]
        data[fname] = np.array(rows)

# rank the keys by best time
sumtimes = [(k,sum(rows.T[-1])) for k, rows in data.items()]
sumtimes.sort(key=lambda x: x[1])
print(sumtimes)

# number of wins
def wins(k):
    num_wins = 0
    for i in range(len(data[k])):
        best_time = data[k][i][-1]
        best_key = k
        for key in data.keys():
            try:
                if data[key][i][-1] < best_time:
                    best_time = data[key][i][-1]
                    best_key = key
            except IndexError:
                pass
        if best_key == k:
            num_wins += 1
    return num_wins

wincount = [(k, wins(k)) for k in data.keys()]
print(wincount)
wincount.sort(key=lambda x: x[1], reverse=True)
