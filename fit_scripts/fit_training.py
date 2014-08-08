import sys
import string

import scipy
from scipy.stats import beta
import numpy


def parse_data(filename):

    x, y = [], []
    f = open(filename, "r")

    for line in f.readlines():
        tokens = string.split(line)
        x.append(float(tokens[0]))
        y.append(float(tokens[1]))

    f.close()

    return x, y

data_file = sys.argv[1]
x, y = parse_data(data_file)

n_bins = 3
bin_step = 100.0/n_bins
bin_edges = numpy.arange(0.0, 100.0 + bin_step, bin_step)

print bin_edges

bin_data = [[] for _ in range(len(bin_edges) - 1)]


for i in range(len(x)):

    n_hb = y[i]
    val  = x[i]

    for j in range(len(bin_data)):

        if (val <= bin_edges[j+1]) and (val > bin_edges[j]):
            bin_data[j].append(n_hb)

for i, data in enumerate(bin_data):

    tdata =  numpy.array(data)

    params = beta.fit(tdata, fscale=1)
    print i, params[0], params[1], params[2], params[3]
    #print tdata





