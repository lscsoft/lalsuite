## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libnoisemodels.la
libnoisemodels_la_SOURCES = \
LALAddVectors.c\
LALColoredNoise.c\
LALCorrelate.c\
LALGEOPsd.c\
LALGaussianNoise.c\
LALGaussianRandomNumber.c\
LALLIGOIPsd.c\
LALNoiseSpectralDensity.c\
LALNormalise.c\
LALInspiralWaveOverlap.c\
LALRandomInspiralSignal.c\
LALTAMAPsd.c\
LALVIRGOPsd.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
