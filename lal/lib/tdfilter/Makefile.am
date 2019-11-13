include $(top_srcdir)/gnuscripts/lalsuite_header_links.am

pkginclude_HEADERS = \
	BandPassTimeSeries.h \
	IIRFilter.h \
	ZPGFilter.h \
	$(END_OF_LIST)

noinst_LTLIBRARIES = libtdfilter.la

libtdfilter_la_SOURCES = \
	BilinearTransform.c \
	CreateZPGFilter.c \
	IIRFilter.c \
	ButterworthTimeSeries.c \
	DestroyIIRFilter.c \
	IIRFilterVector.c \
	CreateIIRFilter.c \
	DestroyZPGFilter.c \
	IIRFilterVectorR.c \
	$(END_OF_LIST)

noinst_HEADERS = \
	ButterworthTimeSeries_source.c \
	CreateIIRFilter_source.c \
	IIRFilterVectorR_source.c \
	IIRFilterVector_source.c \
	$(END_OF_LIST)
