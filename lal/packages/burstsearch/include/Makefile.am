## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = StdBurstSearch.h Thresholds.h ExcessPower.h TFTransform.h BurstSearch.h EPSearch.h Sequence.h TimeSeries.h FrequencySeries.h
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
MAINTAINERCLEANFILES = Sequence.h TimeSeries.h
EXTRA_DIST = SequenceH.m4 Sequence.m4 TimeSeriesH.m4 TimeSeries.m4
all-local: include-link
dvi-local: .dvi-dep
@INCLUDEMKFILE@
@DVIMKFILE@

Sequence.m4 : SequenceH.m4
	touch $@
TimeSeries.m4 : TimeSeriesH.m4
	touch $@
%.h : %.m4
	cd $(srcdir) && m4 $< >$@
