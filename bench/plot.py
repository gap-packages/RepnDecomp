#!/usr/bin/env python3
import numpy as np
import matplotlib, glob

matplotlib.rcParams['pgf.rcfonts'] = False
matplotlib.rcParams['pgf.texsystem'] = 'pdflatex'

import matplotlib.pyplot as plt

NEXT_FIGURE = 1

def do_plot(results, xindex, yindex, xlabel, ylabel, title, fname="graph.png"):
    global NEXT_FIGURE
    to_plot = np.array(results).T
    plt.figure(NEXT_FIGURE)
    NEXT_FIGURE += 1
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.plot(to_plot[xindex], to_plot[yindex], "rx")
    plt.savefig(fname)

# plots a file with cols:
# group size | group id | nr conjugacy classes | degree | time taken
def plot(fname):
    global NEXT_FIGURE
    data = []

    try:
        data = open(fname, "r").read().splitlines()
    except FileNotFoundError:
        return

    pts = np.array([[float(s) for s in line.split(" ")] for line in data])
    to_plot = pts.T

    do_plot(to_plot, 0, 4, "group size", "time taken", f'{fname}: |G| vs t', f'{fname}.png')
    do_plot(to_plot, 2, 4, "num classes", "time taken", f'{fname}: |cc(G)| vs t', f'{fname}.png')
    do_plot(to_plot, 3, 4, "degree", "time taken", f'{fname}: deg vs t', f'{fname}.png')

def tolatex(data, header, output_file):
    output = ""
    output += '\\begin{tabular}{|'
    output += 'n{4}{2}|'*len(data[0])
    output += '}\n'
    output += "\\hline\n"
    output += " & ".join(['{{'+str(col)+'}}' for col in header]) + "\\\\\n"
    output += "\\hline\n"
    output += " \\\\\n".join([" & ".join([str(col) for col in row]) for row in data])
    output += " \\\\\n\\hline"
    output += '\n\\end{tabular}\n'
    open(output_file, "w").write(output)

header = ["|G|", "id", "|cc(G)|", "degree", "time taken"]

for result_csv in glob.glob("*.csv"):
    plot(result_csv)
    with open(result_csv, "r") as result_file:
        lines = result_file.readlines()
        results = [[int(x) for x in row.split(" ")] for row in lines]
        tolatex(results, header, result_csv+".tex")
