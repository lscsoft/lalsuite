## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libnoisemodels.la
libnoisemodels_la_SOURCES = \
LALAddVectors.c\
LALColoredNoise.c\
LALInspiralWaveCorrelate.c\
LALGEOPsd.c\
LALLIGOIPsd.c\
LALNoiseSpectralDensity.c\
LALInspiralComputeChisq.c\
LALInspiralComputeSNRIntegrand.c\
LALInspiralFindEvents.c\
LALInspiralFindEventsCluster.c\
LALInspiralFindLoudestEvent.c\
LALInspiralWaveNormalise.c\
LALInspiralWaveNormaliseLSO.c\
LALInspiralWaveOverlap.c\
LALRandomInspiralSignal.c\
LALStatsREAL4Vector.c\
LALTAMAPsd.c\
LALVIRGOPsd.c\
LALAdvLIGOPsd.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep


@DVIMKFILE@
