## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Units.h DetectorSite.h DetResponse.h Calibration.h \
                     CoherentEstimation.h ResampleTimeSeries.h \
                     LIGOMetadataTables.h LIGOMetadataUtils.h \
                     TemplateBankGeneration.h filters-H1-S2.h  filters-H1-S3.h  \
		     filters-H2-S2.h  filters-H2-S3.h  filters-L1-S2.h  filters-L1-S3.h \
		     Sequence.h TimeSeries.h FrequencySeries.h
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
MAINTAINERCLEANFILES = FrequencySeries.h Sequence.h TimeSeries.h
EXTRA_DIST = FrequencySeries.m4 Sequence.m4 TimeSeries.m4
all-local: include-link
dvi-local: .dvi-dep

FrequencySeries.h : FrequencySeries.m4
	cd $(srcdir) && m4 FrequencySeries.m4 > FrequencySeries.h
Sequence.h : Sequence.m4
	cd $(srcdir) && m4 Sequence.m4 > Sequence.h
TimeSeries.h : TimeSeries.m4
	cd $(srcdir) && m4 TimeSeries.m4 > TimeSeries.h

@INCLUDEMKFILE@
@DVIMKFILE@
