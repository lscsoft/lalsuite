## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libburstsearch.la
libburstsearch_la_SOURCES = SlopeETG.c StdBurstSearch.c TFClustersETG.c AddWhiteNoise.c ComputeExcessPower.c ComputeLikelihood.c ComputeTFPlanes.c CreateTFTiling.c DestroyTFTiling.c PrintTFTileList.c SortTFTiling.c Thresholds.c ComputeFrequencySeries.c CreateComplexDFTParams.c CreateTFPlane.c DestroyComplexDFTParams.c DestroyTFPlane.c FreqSeriesToTFPlane.c TimeSeriesToTFPlane.c EPSearch.c Sequence.c TimeSeries.c FrequencySeries.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
MAINTAINERCLEANFILES = FrequencySeries.c Sequence.c TimeSeries.c
EXTRA_DIST = FrequencySeries.m4 Sequence.m4 TimeSeries.m4
dvi-local: .dvi-dep

FrequencySeries.c: FrequencySeries.m4
	cd $(srcdir) && m4 FrequencySeries.m4 > FrequencySeries.c
Sequence.c: Sequence.m4
	cd $(srcdir) && m4 Sequence.m4 > Sequence.c
TimeSeries.c: TimeSeries.m4
	cd $(srcdir) && m4 TimeSeries.m4 > TimeSeries.c 

@DVIMKFILE@
