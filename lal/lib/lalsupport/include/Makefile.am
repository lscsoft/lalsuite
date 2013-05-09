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
	UserInput.h \
	ConfigFile.h

all-local: include-link
@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
