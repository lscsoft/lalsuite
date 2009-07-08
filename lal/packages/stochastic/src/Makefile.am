## Process this file with automake to produce Makefile.in
lib_LTLIBRARIES = liblalstochastic.la
liblalstochastic_la_SOURCES = StochasticCrossCorrelation.c ZeroPadAndFFT.c \
                           StochasticOptimalFilter.c \
                           StochasticOptimalFilterNormalization.c \
                           OverlapReductionFunction.c StochasticOmegaGW.c \
                           StochasticInverseNoise.c \
                           AstroOmega.c SimulateSB.c \
                           StochasticMC.c SimulatePopcorn.c

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
