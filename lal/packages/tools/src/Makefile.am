## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libtools.la
libtools_la_SOURCES = \
  UnitDefs.c UnitNormalize.c UnitRaise.c UnitMultiply.c UnitCompare.c \
  CreateDetector.c DetResponse.c ComputeTransfer.c CoherentEstimation.c \
  ResampleTimeSeries.c SnglInspiralUtils.c SnglBurstUtils.c \
  BurstUtils.c CoincInspiralUtils.c SimInspiralUtils.c \
  ComputeCalibrationFactors.c LIGOMetadataUtils.c \
  NDTemplateBank.c ComputeStrain.c \
  Sequence.c TimeSeries.c FrequencySeries.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MAINTAINERCLEANFILES = FrequencySeries.c Sequence.c TimeSeries.c
EXTRA_DIST = FrequencySeries.m4 Sequence.m4 TimeSeries.m4 MakeTemplateBank.c
EXTRA_TEXSOURCES = MakeTemplateBank.c
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep

FrequencySeries.c : FrequencySeries.m4
	cd $(srcdir) && @M4@ FrequencySeries.m4 > FrequencySeries.c
Sequence.c : Sequence.m4
	cd $(srcdir) && @M4@ Sequence.m4 > Sequence.c
TimeSeries.c : TimeSeries.m4
	cd $(srcdir) && @M4@ TimeSeries.m4 > TimeSeries.c

@DVIMKFILE@
