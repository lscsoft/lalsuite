import doctest
import functools
import math
import sys
from lal import iterutils
from six.moves import range


#
# check examples in documentation
#


failures = doctest.testmod(iterutils)[0]
if failures:
	sys.exit(bool(failures))


#
# check CDF of randindex.  requires matplotlib for visual check of results.
# without matplotlib, simply runs the code to test for failures
#


def gen_randindex_results(n, N):
	randindex = functools.partial(next, iter(iterutils.randindex(10, 100, n)))
	counts = [0.] * 100
	p = [float("-inf")] * 100
	for i in range(N):
		x, lnP = randindex()
		counts[x] += 1.
		p[x] = lnP
	for i in range(len(counts)):
		counts[i] /= N
		p[i] = math.exp(p[i])
	return list(range(100)), counts, p

try:
	import matplotlib
	matplotlib.use("Agg")
	from matplotlib import figure
	from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
	fig = figure.Figure()
	FigureCanvas(fig)
	axes = fig.gca()
	axes.loglog()
except ImportError:
	axes = None

N = 1000000
x, counts, p = gen_randindex_results(.5, N)
if axes is not None:
	axes.plot(x, counts, "k+", label = "Observed fractions (n=0.5)")
	axes.plot(x, p, "k-", label = "Reported probabilities (n=0.5)")
x, counts, p = gen_randindex_results(1., N)
if axes is not None:
	axes.plot(x, counts, "b+", label = "Observed fractions (n=1)")
	axes.plot(x, p, "b-", label = "Reported probabilities (n=1)")
x, counts, p = gen_randindex_results(2., N)
if axes is not None:
	axes.plot(x, counts, "g+", label = "Observed fractions (n=2)")
	axes.plot(x, p, "g-", label = "Reported probabilities (n=2)")

if axes is not None:
	axes.set_xlim((9., 100.))
	axes.set_title("randindex(10, 100, n), %d samples" % N)
	axes.legend()
	fig.savefig("iterutils_randindex_histogram.png")
