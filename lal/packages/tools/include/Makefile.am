## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = Units.h DetectorSite.h DetResponse.h Calibration.h \
                     CoherentEstimation.h ResampleTimeSeries.h \
                     LIGOMetadataTables.h LIGOMetadataUtils.h \
                     TemplateBankGeneration.h filters-H1-S2.h  filters-H1-S3.h  \
		     filters-H2-S2.h  filters-H2-S3.h  filters-L1-S2.h  filters-L1-S3.h
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
all-local: include-link
dvi-local: .dvi-dep
@INCLUDEMKFILE@
@DVIMKFILE@
