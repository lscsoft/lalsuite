## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libburstsearch.la
libburstsearch_la_SOURCES = PowerETG.c SlopeETG.c StdBurstSearch.c TFClustersETG.c AddWhiteNoise.c ComputeExcessPower.c ComputeLikelihood.c ComputeTFPlanes.c CountEPEvents.c CreateTFTiling.c DestroyTFTiling.c PrintTFTileList.c SortTFTiling.c Thresholds.c ComputeFrequencySeries.c CreateComplexDFTParams.c CreateTFPlane.c DestroyComplexDFTParams.c DestroyTFPlane.c FreqSeriesToTFPlane.c TimeSeriesToTFPlane.c EPSearch.c Sequence.c TimeSeries.c FrequencySeries.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@

Sequence.m4 : SequenceC.m4
	touch $@
TimeSeries.m4 : TimeSeriesC.m4
	touch $@
%.c : %.m4
	cd $(srcdir) && m4 $< >$@
