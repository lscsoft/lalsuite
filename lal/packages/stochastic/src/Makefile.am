## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libstochastic.la
libstochastic_la_SOURCES = StochasticCrossCorrelation.c ZeroPadAndFFT.c \
                           StochasticOptimalFilter.c \
                           StochasticOptimalFilterNormalization.c \
                           OverlapReductionFunction.c StochasticOmegaGW.c \
                           StochasticInverseNoise.c \
                           AstroOmega.c

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
