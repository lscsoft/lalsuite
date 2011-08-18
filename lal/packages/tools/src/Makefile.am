INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

noinst_LTLIBRARIES = libtools.la
libtools_la_SOURCES = \
	CoherentEstimation.c \
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
	LALBarycenter.c \
	LALCalibration.c \
	ResampleTimeSeries.c \
	Segments.c \
	Sequence.c \
	Skymap.c \
	TimeSeries.c \
	UnitCompare.c \
	UnitDefs.c \
	UnitMultiply.c \
	UnitNormalize.c \
	UnitRaise.c

noinst_HEADERS = \
	FrequencySeriesComplex_source.c \
	FrequencySeries_source.c \
	SequenceComplex_source.c \
	Sequence_source.c

MAINTAINERCLEANFILES = \
	TimeSeries.c

EXTRA_DIST = \
	TimeSeries.m4 \
	MakeTemplateBank.c \
	TimeSeriesC.m4

TimeSeries.c : TimeSeries.m4 TimeSeriesC.m4
	cd $(srcdir) && @M4@ TimeSeries.m4 > TimeSeries.c
