#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

NEXT_FIGURE = 1

# plots a file with cols:
# group size | group id | nr conjugacy classes | degree | time taken
def plot(fname):
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
    plt.plot(to_plot[0], to_plot[4], "rx")

    plt.figure(NEXT_FIGURE)
    NEXT_FIGURE += 1
    plt.title("{}: number of classes vs time taken".format(fname))
    plt.plot(to_plot[2], to_plot[4], "rx")

    plt.figure(NEXT_FIGURE)
    NEXT_FIGURE += 1
    plt.title("{}: degree vs time taken".format(fname))
    plt.plot(to_plot[3], to_plot[4], "rx")

if __name__ == "__main__":
    plot("fast.txt")
    plot("serre.txt")
    plt.show()
