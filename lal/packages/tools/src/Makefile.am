INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

noinst_LTLIBRARIES = libtools.la
libtools_la_SOURCES = \
	CoherentEstimation.c \
	CoincInspiralEllipsoid.c \
	CoincInspiralUtils.c \
	CoincRingdownUtils.c \
	ComputeCalibrationFactors.c \
	ComputeDataQualityVector.c \
	ComputeStrain.c \
	ComputeStrainDMT.c \
	ComputeTransfer.c \
	CreateDetector.c \
	DetResponse.c \
	EllipsoidOverlapTools.c \
	FlatLatticeTiling.c \
	FrequencySeries.c \
	GetErrorMatrixFromSnglInspiral.c \
	LALBarycenter.c \
	LALCalibration.c \
	LALTrigScanCluster.c \
	LIGOMetadataUtils.c \
	MultiInspiralUtils.c \
	NDTemplateBank.c \
	ResampleTimeSeries.c \
	Segments.c \
	Sequence.c \
	SimInspiralUtils.c \
	SimRingdownUtils.c \
	Skymap.c \
	SnglBurstUtils.c \
	SnglInspiralUtils.c \
	SnglRingdownUtils.c \
	TimeSeries.c \
	TrigScanEThincaCommon.c \
	UnitCompare.c \
	UnitDefs.c \
	UnitMultiply.c \
	UnitNormalize.c \
	UnitRaise.c

MAINTAINERCLEANFILES = \
	FrequencySeries.c \
	Sequence.c \
	TimeSeries.c
DISTCLEANFILES = *.tex

EXTRA_DIST = \
	FrequencySeries.m4 \
	Sequence.m4 \
	TimeSeries.m4 \
	MakeTemplateBank.c \
	FrequencySeriesC.m4 \
	SequenceC.m4 \
	TimeSeriesC.m4

EXTRA_TEXSOURCES = MakeTemplateBank.c

MOSTLYCLEANFILES = .dvi-dep

dvi-local: .dvi-dep

FrequencySeries.c : FrequencySeries.m4 FrequencySeriesC.m4
	cd $(srcdir) && @M4@ FrequencySeries.m4 > FrequencySeries.c

Sequence.c : Sequence.m4 SequenceC.m4
	cd $(srcdir) && @M4@ Sequence.m4 > Sequence.c

TimeSeries.c : TimeSeries.m4 TimeSeriesC.m4
	cd $(srcdir) && @M4@ TimeSeries.m4 > TimeSeries.c

@DVIMKFILE@
