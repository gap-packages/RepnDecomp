#!/usr/bin/env python3
import numpy as np
import matplotlib, glob, sys, itertools, math

matplotlib.rcParams['pgf.rcfonts'] = False
matplotlib.rcParams['pgf.texsystem'] = 'pdflatex'
matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rcParams['figure.figsize'] = "4,4"

import matplotlib.pyplot as plt

# prettifies a file name for use in a legend
def make_pretty_name(fname):
    s = fname.split('.')[0] # remove .csv
    old_words = s.split('_')[1:] # remove first word
    words = []
    for word in old_words:
        if word == 'mymethod':
            words.append('our method')
        else:
            words.append(word)
    ret = ' '.join([word.title() for word in words])
    return ret

NEXT_FIGURE = 1
log = False

def do_plot(results, xindex, yindex, xlabel, ylabel, title, fname="graph"):
    global NEXT_FIGURE
    to_plot = np.array(results).T
    plt.figure(NEXT_FIGURE)
    NEXT_FIGURE += 1
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.plot(to_plot[xindex], to_plot[yindex], "rx")
    plt.savefig(fname+".pgf")

def do_plot_multiple(results_mult, xindex, yindex, xlabel, ylabel, title, fname, labels):
    global NEXT_FIGURE, log
    marker_color = itertools.cycle(itertools.product("krgbcmy", "x+."))
    plt.figure(NEXT_FIGURE)
    NEXT_FIGURE += 1
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    handles = []
    for i in range(len(results_mult)):
        to_plot = np.array(results_mult[i]).T
        if log:
            for j in range(len(to_plot[yindex])):
                to_plot[yindex][j] = math.log10(to_plot[yindex][j])
        color, marker = next(marker_color)
        handle, = plt.plot(to_plot[xindex], to_plot[yindex], label=labels[i], marker=marker, linestyle="None", color=color)
        handles.append(handle)
    lgd = plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig = plt.gcf()
    fig.savefig(fname+".pgf", bbox_extra_artists=[lgd], bbox_inches='tight')
    fig.savefig(fname+".png", bbox_extra_artists=[lgd], bbox_inches='tight')
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
    do_plot(to_plot, 0, 4, "group size", "time taken", f'{fname}: |G| vs t', f'{fname}_size')
    do_plot(to_plot, 2, 4, "num classes", "time taken", f'{fname}: |cc(G)| vs t', f'{fname}_classes')
    do_plot(to_plot, 3, 4, "degree", "time taken", f'{fname}: deg vs t', f'{fname}_degree')

# plots multiple files on the same graph, for comparison
def plot_multiple(fnames, out_fname):
    time_taken = "time taken (ns)"
    if log:
        time_taken = "log(time taken (ns))"
    ptss = [read_pts(fname) for fname in fnames]
    labels = [make_pretty_name(fname) for fname in fnames]
    do_plot_multiple(ptss, 0, 4, "$|G|$", time_taken, f'', f'{out_fname}_size', labels)
    do_plot_multiple(ptss, 2, 4, "num classes", time_taken, f'', f'{out_fname}_classes', labels)
    do_plot_multiple(ptss, 3, 4, "degree", time_taken, f'', f'{out_fname}_degree', labels)

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

log = sys.argv.pop() == "log"

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
