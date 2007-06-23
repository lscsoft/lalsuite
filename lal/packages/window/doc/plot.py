import math
import matplotlib
matplotlib.rcParams.update({
	"font.size": 8.0,
	"axes.titlesize": 10.0,
	"axes.labelsize": 10.0,
	"xtick.labelsize": 8.0,
	"ytick.labelsize": 8.0,
	"legend.fontsize": 8.0,
	"figure.dpi": 300,
	"savefig.dpi": 300,
	"text.usetex": True     # render all text with TeX
})
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import numpy

from pylal.window import *

def make_windows(n, kaiser_beta, creighton_beta, tukey_beta, gauss_beta):
	return {
		"rectangle": XLALCreateRectangularREAL8Window(n),
		"Hann": XLALCreateHannREAL8Window(n),
		"Welch": XLALCreateWelchREAL8Window(n),
		"Bartlett": XLALCreateBartlettREAL8Window(n),
		"Parzen": XLALCreateParzenREAL8Window(n),
		"Papoulis": XLALCreatePapoulisREAL8Window(n),
		"Hamming": XLALCreateHammingREAL8Window(n),
		"Kaiser": XLALCreateKaiserREAL8Window(n, kaiser_beta),
		"Creighton": XLALCreateCreightonREAL8Window(n, creighton_beta),
		"Tukey": XLALCreateTukeyREAL8Window(n, tukey_beta),
		"Gauss": XLALCreateGaussREAL8Window(n, gauss_beta)
	}

x = numpy.arange(1001, dtype = "float64") / 500.5 - 1

fig = figure.Figure()
FigureCanvas(fig)
fig.set_size_inches(6.5, 6.5 / ((1 + math.sqrt(5)) / 2))
axes = fig.gca()
axes.grid(True)
axes.set_xlabel("$y$")
axes.set_ylabel("$w(y)$")

windows = make_windows(1001, 6, 2, 0.5, 3)

for window in windows.values():
	axes.plot(x, window.data)
axes.legend(windows.keys())

fig.savefig("window_t.eps")
