#!/usr/bin/env sage

# trying to load() this script would not work at Sage prompt
load("crossing.sage")

if __name__ == "__main__":
    import sys
    m = int(sys.argv[1])
    compute_alpha(m)
