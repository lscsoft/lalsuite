## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = StdBurstSearch.h Thresholds.h ExcessPower.h TFTransform.h EPSearch.h Sequence.h TimeSeries.h FrequencySeries.h
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
MAINTAINERCLEANFILES = FrequencySeries.h Sequence.h TimeSeries.h
EXTRA_DIST = FrequencySeriesH.m4 FrequencySeries.m4 SequenceH.m4 Sequence.m4 TimeSeriesH.m4 TimeSeries.m4
all-local: include-link
dvi-local: .dvi-dep
@INCLUDEMKFILE@
@DVIMKFILE@

FrequencySeries.m4 : FrequencySeriesH.m4
	cd $(srcdir) && touch $@
Sequence.m4 : SequenceH.m4
	cd $(srcdir) && touch $@
TimeSeries.m4 : TimeSeriesH.m4
	cd $(srcdir) && touch $@
%.h : %.m4
	cd $(srcdir) && m4 $< >$@
