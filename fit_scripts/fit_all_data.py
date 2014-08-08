import sys
import string
import math
import os.path
import cPickle
import Bio.PDB
import numpy
import scipy
from scipy.stats import beta as beta_dist
from matplotlib import pyplot
from scipy.special import binom, gamma
from matplotlib import pyplot

def load_pickle(filename):
    f = open(filename,"rb")
    p = cPickle.load(f)
    f.close()
    return(p)



def parse_data(filename):

    data = []

    f = open(filename, "r")
    lines = f.readlines()

    for line in lines:
        if line[0] == "#" or len(line) < 4:
            continue

        tokens = string.split(line)

        start = int(tokens[0])
        end = int(tokens[1])
        uptake = float(tokens[2])/100.0

        residues = range(start, end + 1)
        length = len(residues)

        data_point = dict()

        data_point['residues'] = residues
        data_point['uptake'] = uptake
        data_point['length'] = length

        data.append(data_point)

    f.close()

    return data



def factorial(x):

    try:
        return gamma(x + 1)
    except OverflowError:
        print "Overflow, x =",x
        exit(0)



def B(x, y):
    return factorial(x - 1) * factorial(y - 1) / factorial(x + y - 1)



def get_uptake_avg_hbonds(pickle_filenames, data):


    uptake = []
    avg_hbonds = []

    for pickle_filename in pickle_filenames:

        hbonds = load_pickle(pickle_filename)



        for d in data:

            #print d
            n_bonds = 0.0

            for residue in d['residues']:
                if residue in hbonds:
                    n_bonds += 1

            # print d, d['uptake'], n_bonds/d['length'], d['length']

            if n_bonds/d['length'] < 0.1 and d['uptake'] < 0.60:
                print "Outlier:", d, d['uptake'], n_bonds/d['length'], pickle_filename
                continue

            uptake.append(d['uptake'] * 100.0)
            avg_hbonds.append(n_bonds/d['length'])

    return uptake, avg_hbonds



def beta_binom(k, n, aplha, beta):
    return binom(n, k) * B(k + alpha, n - k + beta) / B(alpha, beta)


#for k in range(0, n + 1):
#   print "P(N =%3i) = %6.4f" % (k, beta_binom(k))






if __name__ == "__main__":

    data_filename = sys.argv[1]
    data = parse_data(data_filename)

    x = []
    y = []



    pickle_filenames = ["clean_structures/1axi.pqr.cpickle",
                        "clean_structures/3hhr.pqr.cpickle",
                        "clean_structures/1hgu.pqr.cpickle",
                        "clean_structures/hgh_model.cpickle"]


    uptakes, avg_hbonds = get_uptake_avg_hbonds(pickle_filenames, data)

    strand_length = 5
    data_bins = 4
    bin_edges = numpy.linspace(0.0, 1.0, num=data_bins + 1, endpoint=True)
    print bin_edges
    bin_data = [[] for _ in range(data_bins)]

    for i, uptake in enumerate(uptakes):

        n_hb = avg_hbonds[i]
        uptake /= 100.0

        for j in range(data_bins):

            if (uptake <= bin_edges[j+1]) and (uptake > bin_edges[j]):
                bin_data[j].append(n_hb)

    numpy.set_printoptions(linewidth=1000000)



    probs = []
    pvals = []

    pyplot.figure(figsize=(15,10))

    for i, data in enumerate(bin_data):

        print i
        tdata =  numpy.array(data)

        params = beta_dist.fit(tdata, floc=0) 

        (alpha, beta, floc, fshape) = params

        print alpha, beta

        vals = numpy.arange(0, strand_length + 1)
        fit_hist = numpy.array([beta_binom(val, strand_length, alpha, beta) for val in vals])

        #data_bin_edges = numpy.linspace(0.0, 1.0, num=strand_length + 1, endpoint=True)

        data_bin_edges = numpy.linspace(-0.5, strand_length + 0.5, num=strand_length + 2,  endpoint=True) / float(strand_length)

        #print data_bin_edges


        data_hist, data_bins = numpy.histogram(tdata, data_bin_edges)

        data_hist = numpy.array(data_hist) / float(sum(data_hist))

        #print len(fit_hist), fit_hist
        #print len(data_hist), data_hist

        pyplot.subplot(2,math.ceil(len(bin_data)/2.0), i+1)

        bar_width = 0.35

        pyplot.title("Uptake %i - %i, (a = %1.2f, b = %1.2f)" % (int(bin_edges[i] * 100), 
                            int(bin_edges[i+1] * 100), alpha, beta))
        pyplot.bar(vals + bar_width/2, fit_hist, bar_width, color = 'DarkSlateBlue', alpha=0.6)
        pyplot.bar(vals + bar_width/2+bar_width, data_hist, bar_width, color = "r", alpha=0.6)
        pyplot.xticks(vals + 1.5*bar_width, vals, fontsize=16)
        pyplot.ylim([0, 1.0])

    pyplot.savefig("test_fit.png")
