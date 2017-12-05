include $(top_srcdir)/gnuscripts/lalsuite_header_links.am

pkginclude_HEADERS = \
	Calibration.h \
	CoherentEstimation.h \
	ComputeDataQualityVector.h \
	DetResponse.h \
	DetectorSite.h \
	EllipsoidOverlapTools.h \
	FrequencySeries.h \
	LALCosmologyCalculator.h \
	LALDetectors.h \
	LALDict.h \
	LALList.h \
	LALValue.h \
	ResampleTimeSeries.h \
	Segments.h \
	Sequence.h \
	Skymap.h \
	TimeSeries.h \
	TimeSeriesInterp.h \
	TriggerInterpolation.h \
	Units.h \
	$(END_OF_LIST)

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
	FrequencySeries.c \
	LALCosmologyCalculator.c \
	LALDict.c \
	LALList.c \
	LALValue.c \
	ResampleTimeSeries.c \
	Segments.c \
	Sequence.c \
	Skymap.c \
	TimeSeries.c \
	TimeSeriesInterp.c \
	TriggerInterpolation.c \
	UnitCompare.c \
	UnitDefs.c \
	UnitMultiply.c \
	UnitNormalize.c \
	UnitRaise.c \
	$(END_OF_LIST)

noinst_HEADERS = \
	FrequencySeriesComplex_source.c \
	FrequencySeries_source.c \
	LALValue_private.h \
	SequenceComplex_source.c \
	Sequence_source.c \
	TimeSeries_source.c \
	$(END_OF_LIST)

EXTRA_DIST = \
	MakeTemplateBank.c \
	$(END_OF_LIST)
