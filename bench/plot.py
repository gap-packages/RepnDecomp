#!/usr/bin/env python3
import numpy as np
import matplotlib, glob, sys, itertools

matplotlib.rcParams['pgf.rcfonts'] = False
matplotlib.rcParams['pgf.texsystem'] = 'pdflatex'
matplotlib.rcParams['figure.figsize'] = 10,10

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

marker_color = itertools.cycle(itertools.product("xo+", "rgbcmyk"))

def do_plot_multiple(results_mult, xindex, yindex, xlabel, ylabel, title, fname, labels):
    global NEXT_FIGURE
    plt.figure(NEXT_FIGURE)
    NEXT_FIGURE += 1
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    handles = []
    for i in range(len(results_mult)):
        to_plot = np.array(results_mult[i]).T
        marker, color = next(marker_color)
        handle, = plt.plot(to_plot[xindex], to_plot[yindex], label=labels[i], marker=marker, linestyle="None", color=color)
        handles.append(handle)
    lgd = plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig(fname, bbox_extra_artists=[lgd], bbox_inches='tight')
    print(fname)

def read_pts(fname):
    data = []

    try:
        data = open(fname, "r").read().splitlines()
    except FileNotFoundError:
        return

    return np.array([[float(s) for s in line.split(" ")] for line in data])

# plots a file with cols:
# group size | group id | nr conjugacy classes | degree | time taken
def plot(fname):
    pts = read_pts(fname)
    to_plot = pts.T
    do_plot(to_plot, 0, 4, "group size", "time taken", f'{fname}: |G| vs t', f'{fname}_size.png')
    do_plot(to_plot, 2, 4, "num classes", "time taken", f'{fname}: |cc(G)| vs t', f'{fname}_classes.png')
    do_plot(to_plot, 3, 4, "degree", "time taken", f'{fname}: deg vs t', f'{fname}_degree.png')

# plots multiple files on the same graph, for comparison
def plot_multiple(fnames, out_fname):
    ptss = [read_pts(fname) for fname in fnames]
    do_plot_multiple(ptss, 0, 4, "group size", "time taken", f'{out_fname}: |G| vs t', f'{out_fname}_size.png', fnames)
    do_plot_multiple(ptss, 2, 4, "num classes", "time taken", f'{out_fname}: |cc(G)| vs t', f'{out_fname}_classes.png', fnames)
    do_plot_multiple(ptss, 3, 4, "degree", "time taken", f'{out_fname}: deg vs t', f'{out_fname}_degree.png', fnames)

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

fnames = sys.argv[2:]

if len(fnames) > 1:
    plot_multiple(fnames, sys.argv[1])
else:
    plot(fnames[0])

for result_csv in fnames:
    with open(result_csv, "r") as result_file:
        lines = result_file.readlines()
        results = [[int(x) for x in row.split(" ")] for row in lines]
        tolatex(results, header, result_csv+".tex")
