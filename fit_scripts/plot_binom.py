import numpy
from scipy.special import binom, gamma

from matplotlib import pyplot
import sys

def factorial(x):

    try:
        return gamma(x + 1)
    except OverflowError:
        print "Overflow, x =",x
        exit(0)

def B(x, y):
    return factorial(x - 1) * factorial(y - 1) / factorial(x + y - 1)


n      = int(sys.argv[1])
mu     = float(sys.argv[2])
sigma  = float(sys.argv[3])

alpha = - mu * (mu * mu - mu * n + sigma * sigma) / (sigma * sigma * n + mu * mu - mu * n)
beta  = (n * alpha) / mu - alpha


alpha = float(sys.argv[2])
beta  = float(sys.argv[3])

if (alpha < 0.0) or (beta < 0.0):
    print "ERROR: Negative parameter value:"
    print "alpha =", alpha, "beta =", beta
    exit(0)

sigma = numpy.sqrt( n * alpha * beta * (alpha + beta + n) / ((alpha + beta) * (alpha + beta) * (1 + alpha + beta)))
mu = n * alpha / (alpha + beta)

print "alpha =", alpha, "beta =", beta
print "mu = %f sigma = %f" % (mu, sigma)


def beta_binom(k):
    return binom(n, k) * B(k + alpha, n - k + beta) / B(alpha, beta)


for k in range(0, n + 1):
    print "P(N =%3i) = %6.4f" % (k, beta_binom(k))


pyplot.rc('text', usetex=True)
pyplot.rc('font', family='serif')


vals = numpy.arange(0, n + 1)
probs = numpy.array([beta_binom(val) for val in vals])

bar_width = 0.55
pyplot.bar(vals + bar_width/2, probs, bar_width, color = 'DarkSlateBlue', alpha=0.6)

pyplot.title(r"$n = %i,\ \mu= %5.2f,\ \sigma = %5.2f\ (\alpha = %5.2f,\ \beta = %5.2f)$" % (n, mu, sigma, alpha, beta), fontsize=20)


val_texts = [r"$%i$" % (val) for val in vals]
pyplot.xlabel(r"$k$", fontsize=16)
pyplot.xticks(vals + bar_width, val_texts, fontsize=16)
pyplot.xlim([0.0, numpy.amax(vals) + bar_width*2])

pyplot.yticks(fontsize=16)
pyplot.ylabel(r"$P(N_\mathrm{HB}=k)$", fontsize=16)

pyplot.grid(True)
pyplot.savefig("bar.png")
