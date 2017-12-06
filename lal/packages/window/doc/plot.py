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
from scipy import fftpack

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

def plot_windows_t(x, windows):
	fig = figure.Figure()
	FigureCanvas(fig)
	fig.set_size_inches(6.5, 6.5 / ((1 + math.sqrt(5)) / 2))
	axes = fig.gca()
	axes.grid(True)
	axes.set_xlabel(r"$y$")
	axes.set_ylabel(r"$w(y)$")

	for window in windows.values():
		axes.plot(x, window.data)
	axes.legend(windows.keys())

	return fig

def plot_windows_f(x, windows):
	fig = figure.Figure()
	FigureCanvas(fig)
	fig.set_size_inches(6.5, 6.5 / ((1 + math.sqrt(5)) / 2))
	axes = fig.gca()
	axes.grid(True)
	axes.loglog()
	axes.set_xlabel(r"$y^{-1}$")
	axes.set_ylabel(r"$\left| \tilde{w}(y^{-1}) \right| / L$")

	L = len(x)
	pad_factor = 100
	padding = numpy.zeros(((pad_factor - 1) * L,))
	for window in windows.values():
		# get window data, and pad with zeros to make its length
		# pad_factor * L
		w = numpy.hstack((window.data, padding))

		# Fourier transform, and extract the bit to plot
		w = fftpack.fft(w)[:4 * L]

		# compute the horizontal co-ordinate
		yinv = numpy.arange(0, len(w), dtype = "float64") / pad_factor

		# plot
		axes.plot(yinv, numpy.abs(w) / L)
	axes.set_xlim([1e-1, 1e+2])
	axes.set_ylim([1e-7, 1e+0])
	axes.legend(windows.keys())

	return fig

L = 1001
x = numpy.arange(L) / ((L - 1) / 2.0) - 1
windows = make_windows(L, 6, 2, 0.5, 3)

plot_windows_t(x, windows).savefig("window_t.eps")
plot_windows_f(x, windows).savefig("window_f.eps")
