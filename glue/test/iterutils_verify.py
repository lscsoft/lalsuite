import doctest
from glue import iterutils

doctest.testmod(iterutils)


#
# check CDF of randindex
#


import math
import matplotlib
matplotlib.use("Agg")
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

def gen_randindex_results(n, N):
	randindex = iter(iterutils.randindex(10, 100, n)).next
	counts = [0.] * 100
	p = [float("-inf")] * 100
	for i in xrange(N):
		x, lnP = randindex()
		counts[x] += 1.
		p[x] = lnP
	for i in xrange(len(counts)):
		counts[i] /= N
		p[i] = math.exp(p[i])
	return range(100), counts, p

fig = figure.Figure()
FigureCanvas(fig)
axes = fig.gca()
axes.loglog()

N = 1000000
x, counts, p = gen_randindex_results(.5, N)
axes.plot(x, counts, "k+", label = "Observed fractions (n=0.5)")
axes.plot(x, p, "k-", label = "Reported probabilities (n=0.5)")
x, counts, p = gen_randindex_results(1., N)
axes.plot(x, counts, "b+", label = "Observed fractions (n=1)")
axes.plot(x, p, "b-", label = "Reported probabilities (n=1)")
x, counts, p = gen_randindex_results(2., N)
axes.plot(x, counts, "g+", label = "Observed fractions (n=2)")
axes.plot(x, p, "g-", label = "Reported probabilities (n=2)")

axes.set_xlim((9., 100.))
axes.set_title("randindex(10, 100), %d samples" % N)
axes.legend()
fig.savefig("iterutils_randindex_histogram.png")
