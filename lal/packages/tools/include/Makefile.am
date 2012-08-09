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

MAINTAINERCLEANFILES =

all-local: include-link

@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
