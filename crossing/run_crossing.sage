#!/usr/bin/env sage

# trying to load() this script would not work at Sage prompt
load("crossing.sage")

if __name__ == "__main__":
    import sys, sage.interfaces.gap

    # needed for alpha_7 and higher
    sage.interfaces.gap.set_gap_memory_pool_size(10000000000)

    m = int(sys.argv[1])
    compute_alpha(m)
