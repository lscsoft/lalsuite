pkginclude_HEADERS = \
	Audio.h \
	FileIO.h \
	LALCache.h \
	LALInitBarycenter.h \
	LALMathematica.h \
	LogPrintf.h \
	PrintFTSeries.h \
	PrintVector.h \
	ReadFTSeries.h \
	ReadFiltersFile.h \
	ReadNoiseSpectrum.h \
	SegmentsIO.h \
	StreamInput.h \
	StreamOutput.h \
	UserInput.h ConfigFile.h

CLEANFILES =
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
all-local: include-link
dvi-local: .dvi-dep
@INCLUDEMKFILE@
@DVIMKFILE@

swig_headers = $(pkginclude_HEADERS)
CLEANFILES += $(swig_clean)
@SWIG_HEADER_MK@
