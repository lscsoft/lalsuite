## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libtools.la
libtools_la_SOURCES = \
  UnitDefs.c UnitNormalize.c UnitRaise.c UnitMultiply.c UnitCompare.c \
  CreateDetector.c DetResponse.c LALCalibration.c ComputeTransfer.c \
  CoherentEstimation.c ResampleTimeSeries.c SnglInspiralUtils.c \
  MultiInspiralUtils.c SnglBurstUtils.c SnglRingdownUtils.c \
  CoincInspiralUtils.c CoincRingdownUtils.c SimInspiralUtils.c \
  ComputeCalibrationFactors.c LIGOMetadataUtils.c \
  NDTemplateBank.c ComputeStrain.c ComputeStrainDMT.c \
  Sequence.c TimeSeries.c FrequencySeries.c \
  Segments.c LALTrigScanCluster.c TrigScanEThincaCommon.c \
  EllipsoidOverlapTools.c GetErrorMatrixFromSnglInspiral.c \
  CoincInspiralEllipsoid.c SimRingdownUtils.c \
  FlatLatticeTiling.c ComputeDataQualityVector.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MAINTAINERCLEANFILES = FrequencySeries.c Sequence.c TimeSeries.c
EXTRA_DIST = \
  FrequencySeries.m4 Sequence.m4 TimeSeries.m4 MakeTemplateBank.c \
  FrequencySeriesC.m4 SequenceC.m4 TimeSeriesC.m4
EXTRA_TEXSOURCES = MakeTemplateBank.c
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep

FrequencySeries.c : FrequencySeries.m4 FrequencySeriesC.m4
	cd $(srcdir) && @M4@ FrequencySeries.m4 > FrequencySeries.c
Sequence.c : Sequence.m4 SequenceC.m4
	cd $(srcdir) && @M4@ Sequence.m4 > Sequence.c
TimeSeries.c : TimeSeries.m4 TimeSeriesC.m4
	cd $(srcdir) && @M4@ TimeSeries.m4 > TimeSeries.c

@DVIMKFILE@
