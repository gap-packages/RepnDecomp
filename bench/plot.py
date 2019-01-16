#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

NEXT_FIGURE = 1

def plot(fname):
    global NEXT_FIGURE
    data = open(fname, "r").read().splitlines()
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

if __name__ == "__main__":
    plot("fast.txt")
    plot("serre.txt")
    plt.show()
