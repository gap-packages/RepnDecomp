#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

NEXT_FIGURE = 1

# plots a file with cols "group size | nr classes | time taken"
def plot_size(fname):
    global NEXT_FIGURE
    data = []

    try:
        data = open(fname, "r").read().splitlines()
    except FileNotFoundError:
        return

    pts = np.array([[int(s) for s in line.split(" ")] for line in data])
    to_plot = pts.T

    plt.figure(NEXT_FIGURE)
    NEXT_FIGURE += 1
    plt.title("{}: group size vs time taken".format(fname))
    plt.plot(to_plot[0], to_plot[2], "rx")

    plt.figure(NEXT_FIGURE)
    NEXT_FIGURE += 1
    plt.title("{}: number of classes vs time taken".format(fname))
    plt.plot(to_plot[1], to_plot[2], "rx")

# plots a file with cols "degree | time taken". group size is assumed
# fixed
def plot_degree(fname):
    global NEXT_FIGURE
    data = []

    try:
        data = open(fname, "r").read().splitlines()
    except FileNotFoundError:
        return

    pts = np.array([[int(s) for s in line.split(" ")] for line in data])
    to_plot = pts.T

    plt.figure(NEXT_FIGURE)
    NEXT_FIGURE += 1
    plt.title("{}: degree vs time taken".format(fname))
    plt.plot(to_plot[0], to_plot[1], "rx")

if __name__ == "__main__":
    plot_size("fast.txt")
    plot_size("serre.txt")
    plot_size("canonical_serre.txt")
    plot_size("canonical_fast.txt")
    plot_degree("fast_deg.txt")
    plot_degree("serre_deg.txt")
    plt.show()
