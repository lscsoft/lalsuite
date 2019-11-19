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
	ReadNoiseSpectrum.h \
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
	ReadFrequencySeries.c \
	ReadTimeSeries.c \
	ReadNoiseSpectrum.c \
	StreamSeriesInput.c \
	StreamSeriesOutput.c \
	StreamSequenceInput.c \
	StreamVectorInput.c \
	StreamVectorSequenceInput.c \
	TranslateAngles.c \
	TranslateMJD.c \
	UserInput.c \
	UserInputParse.c \
	UserInputPrint.c \
	LALgetopt.c \
	$(END_OF_LIST)

liblalsupport_la_LIBADD = \
	libtypeZ.la \
	libtypeC.la \
	libtypeD.la \
	libtypeS.la \
	libtypeI2.la \
	libtypeI4.la \
	libtypeI8.la \
	libtypeU2.la \
	libtypeU4.la \
	libtypeU8.la \
	libtypeCHAR.la \
	libtype.la

noinst_LTLIBRARIES = $(liblalsupport_la_LIBADD)

per_type_sources = PrintVector.c

libtypeZ_la_SOURCES = $(per_type_sources)
libtypeZ_la_CPPFLAGS = $(AM_CPPFLAGS) \
	-DTYPECODE='Z' \
	-DTYPE='COMPLEX16' \
	-DFMT='"%i %g %g\n"' -DARG='creal(vector->data[i]),cimag(vector->data[i])'

libtypeC_la_SOURCES = $(per_type_sources)
libtypeC_la_CPPFLAGS = $(AM_CPPFLAGS) \
	-DTYPECODE='C' \
	-DTYPE='COMPLEX8' \
	-DFMT='"%i %g %g\n"' -DARG='crealf(vector->data[i]),cimagf(vector->data[i])'

libtypeD_la_SOURCES = $(per_type_sources)
libtypeD_la_CPPFLAGS = $(AM_CPPFLAGS) \
	-DTYPECODE='D' \
	-DTYPE='REAL8' \
	-DFMT='"%i %g\n"' -DARG='vector->data[i]'

libtypeS_la_SOURCES = $(per_type_sources)
libtypeS_la_CPPFLAGS = $(AM_CPPFLAGS) \
	-DTYPECODE='S' \
	-DTYPE='REAL4' \
	-DFMT='"%i %g\n"' -DARG='vector->data[i]'

libtypeI2_la_SOURCES = $(per_type_sources)
libtypeI2_la_CPPFLAGS = $(AM_CPPFLAGS) \
	-DTYPECODE='I2' \
	-DTYPE='INT2' \
	-DFMT='"%i %i\n"' -DARG='vector->data[i]'

libtypeI4_la_SOURCES = $(per_type_sources)
libtypeI4_la_CPPFLAGS = $(AM_CPPFLAGS) \
	-DTYPECODE='I4' \
	-DTYPE='INT4' \
	-DFMT='"%i %i\n"' -DARG='vector->data[i]'

libtypeI8_la_SOURCES = $(per_type_sources)
libtypeI8_la_CPPFLAGS = $(AM_CPPFLAGS) \
	-DTYPECODE='I8' \
	-DTYPE='INT8' \
	-DFMT='"%i %0.0f\n"' -DARG='(REAL8)vector->data[i]'

libtypeU2_la_SOURCES = $(per_type_sources)
libtypeU2_la_CPPFLAGS = $(AM_CPPFLAGS) \
	-DTYPECODE='U2' \
	-DTYPE='UINT2' \
	-DFMT='"%i %i\n"' -DARG='vector->data[i]'

libtypeU4_la_SOURCES = $(per_type_sources)
libtypeU4_la_CPPFLAGS = $(AM_CPPFLAGS) \
	-DTYPECODE='U4' \
	-DTYPE='UINT4' \
	-DFMT='"%i %i\n"' -DARG='vector->data[i]'

libtypeU8_la_SOURCES = $(per_type_sources)
libtypeU8_la_CPPFLAGS = $(AM_CPPFLAGS) \
	-DTYPECODE='U8' \
	-DTYPE='UINT8' \
	-DFMT='"%i %0.0f\n"' -DARG='(REAL8)vector->data[i]'

libtypeCHAR_la_SOURCES = $(per_type_sources)
libtypeCHAR_la_CPPFLAGS = $(AM_CPPFLAGS) \
	-DTYPECODE='CHAR' \
	-DTYPE='CHAR' \
	-DFMT='"%i %c\n"' -DARG='vector->data[i]'

libtype_la_SOURCES = $(per_type_sources)
libtype_la_CPPFLAGS = $(AM_CPPFLAGS) \
	-DTYPECODE='' \
	-DTYPE='REAL4' \
	-DFMT='"%i %f\n"' -DARG='vector->data[i]'

liblalsupport_la_CFLAGS = $(AM_CFLAGS) $(ZLIB_CFLAGS) $(HDF5_CFLAGS)
liblalsupport_la_CPPFLAGS = $(AM_CPPFLAGS) $(HDF5_CPPFLAGS) $(PAGER_CPPFLAGS)
liblalsupport_la_LDFLAGS = ../liblal.la $(AM_LDFLAGS) $(HDF5_LDFLAGS) $(ZLIB_LIBS) $(HDF5_LIBS) -version-info $(LIBVERSION_SUPPORT)

noinst_HEADERS = \
	H5FileIOArrayHL_source.c \
	H5FileIOArray_source.c \
	H5FileIOFrequencySeries_source.c \
	H5FileIOScalar_source.c \
	H5FileIOTimeSeries_source.c \
	H5FileIOVectorHL_source.c \
	H5FileIOVector_source.c \
	H5ThreadSafe.c \
	PrintFrequencySeries_source.c \
	PrintTimeSeries_source.c \
	ReadFrequencySeries_source.c \
	ReadTimeSeries_source.c \
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
