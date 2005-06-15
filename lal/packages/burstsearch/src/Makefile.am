## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libburstsearch.la
libburstsearch_la_SOURCES = SlopeETG.c StdBurstSearch.c TFClustersETG.c AddWhiteNoise.c ComputeExcessPower.c ComputeLikelihood.c ComputeTFPlanes.c CreateTFTiling.c DestroyTFTiling.c PrintTFTileList.c SortTFTiling.c Thresholds.c ComputeFrequencySeries.c CreateTFPlane.c DestroyTFPlane.c FreqSeriesToTFPlane.c EPSearch.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
