pkginclude_HEADERS = \
	Calibration.h \
	CoherentEstimation.h \
	ComputeDataQualityVector.h \
	DetResponse.h \
	DetectorSite.h \
	EllipsoidOverlapTools.h \
	FlatLatticeTiling.h \
	FrequencySeries.h \
	LALBarycenter.h \
	LALCalibration.h \
	LALDetectors.h \
	ResampleTimeSeries.h \
	Segments.h \
	Sequence.h \
	Skymap.h \
	TimeSeries.h \
	Units.h

CLEANFILES =

MAINTAINERCLEANFILES = \
	TimeSeries.h

EXTRA_DIST = \
	TimeSeries.m4 \
	TimeSeriesH.m4

all-local: include-link

TimeSeries.h : TimeSeries.m4 TimeSeriesH.m4
	cd $(srcdir) && @M4@ TimeSeries.m4 > TimeSeries.h

@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
include $(top_srcdir)/swig/swig-header.am
