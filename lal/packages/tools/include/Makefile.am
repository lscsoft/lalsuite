## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Units.h LALDetectors.h DetectorSite.h DetResponse.h \
		     LALCalibration.h Calibration.h CoherentEstimation.h \
		     ResampleTimeSeries.h LIGOMetadataTables.h \
		     LIGOMetadataUtils.h TemplateBankGeneration.h \
		     Sequence.h TimeSeries.h FrequencySeries.h \
		     Segments.h LALTrigScanCluster.h EllipsoidOverlapTools.h \
		     CoincInspiralEllipsoid.h FlatLatticeTiling.h \
		     ComputeDataQualityVector.h TrigScanEThincaCommon.h
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
MAINTAINERCLEANFILES = FrequencySeries.h Sequence.h TimeSeries.h
EXTRA_DIST = \
  FrequencySeries.m4 Sequence.m4 TimeSeries.m4 \
  FrequencySeriesH.m4 SequenceH.m4 TimeSeriesH.m4
all-local: include-link
dvi-local: .dvi-dep

FrequencySeries.h : FrequencySeries.m4 FrequencySeriesH.m4
	cd $(srcdir) && @M4@ FrequencySeries.m4 > FrequencySeries.h
Sequence.h : Sequence.m4 SequenceH.m4
	cd $(srcdir) && @M4@ Sequence.m4 > Sequence.h
TimeSeries.h : TimeSeries.m4 TimeSeriesH.m4
	cd $(srcdir) && @M4@ TimeSeries.m4 > TimeSeries.h

@INCLUDEMKFILE@
@DVIMKFILE@
