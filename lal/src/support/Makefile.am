include $(top_srcdir)/gnuscripts/lalsuite_header_links.am

pkginclude_HEADERS = \
	Audio.h \
	FileIO.h \
	H5FileIO.h \
	LALCache.h \
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
	TranslateAngles.h \
	TranslateMJD.h \
	UserInput.h \
	UserInputParse.h \
	UserInputPrint.h \
	ConfigFile.h \
	LALgetopt.h \
	$(END_OF_LIST)

lib_LTLIBRARIES = liblalsupport.la

liblalsupport_la_SOURCES = \
	Audio.c \
	ConfigFile.c \
	FileIO.c \
	H5FileIOHighLevel.c \
	H5FileIOLowLevel.c \
	H5FileIOMidLevel.c \
	LALCache.c \
	LALMath3DPlot.c \
	LALMathNDPlot.c \
	LogPrintf.c \
	PrintFrequencySeries.c \
	PrintTimeSeries.c \
	PrintVector.c \
	ReadFiltersFile.c \
	ReadFrequencySeries.c \
	ReadTimeSeries.c \
	ReadNoiseSpectrum.c \
	SegmentsIO.c \
	StreamSeriesInput.c \
	StreamSeriesOutput.c \
	StreamSequenceInput.c \
	StreamGridInput.c \
	StreamGridOutput.c \
	StreamVectorInput.c \
	StreamVectorSequenceInput.c \
	TranslateAngles.c \
	TranslateMJD.c \
	UserInput.c \
	UserInputParse.c \
	UserInputPrint.c \
	LALgetopt.c \
	$(END_OF_LIST)

liblalsupport_la_LIBADD = ../liblal.la $(ZLIB_LIBS) $(HDF5_LIBS)
liblalsupport_la_CFLAGS = $(AM_CFLAGS) $(ZLIB_CFLAGS) $(HDF5_CFLAGS)
liblalsupport_la_CPPFLAGS = $(AM_CPPFLAGS) $(HDF5_CPPFLAGS) $(PAGER_CPPFLAGS)
liblalsupport_la_LDFLAGS = $(AM_LDFLAGS) $(HDF5_LDFLAGS) -version-info $(LIBVERSION_SUPPORT)

noinst_HEADERS = \
	H5FileIOArrayHL_source.c \
	H5FileIOArray_source.c \
	H5FileIOFrequencySeries_source.c \
	H5FileIOScalar_source.c \
	H5FileIOTimeSeries_source.c \
	H5FileIOVectorHL_source.c \
	H5FileIOVector_source.c \
	PrintFrequencySeries_source.c \
	PrintTimeSeries_source.c \
	PrintVector_source.c \
	ReadFrequencySeries_source.c \
	ReadTimeSeries_source.c \
	StreamGridInput_source.c \
	StreamGridOutput_source.c \
	StreamSequenceInput_source.c \
	StreamSequenceInputComplex_source.c \
	StreamSeriesInput_source.c \
	StreamSeriesOutput_source.c \
	StreamVectorInput_source.c \
	StreamVectorSequenceInput_source.c \
	$(END_OF_LIST)

EXTRA_DIST = \
	LALMathematicaHplot1.png \
	$(END_OF_LIST)
