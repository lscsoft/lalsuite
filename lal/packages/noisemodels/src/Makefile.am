## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libnoisemodels.la
libnoisemodels_la_SOURCES = \
LALAddVectors.c\
LALColoredNoise.c\
LALGEOPsd.c\
LALEGOPsd.c\
LALLIGOIPsd.c\
LALNoiseSpectralDensity.c\
LALStatsREAL4Vector.c\
LALTAMAPsd.c\
LALVIRGOPsd.c\
LALAdvLIGOPsd.c
INCLUDES = -I$(top_builddir)/include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep


@DVIMKFILE@
