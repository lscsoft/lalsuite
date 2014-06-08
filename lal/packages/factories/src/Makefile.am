include $(top_srcdir)/gnuscripts/lallibs_header_links.am

AM_CPPFLAGS = -I$(top_builddir)/include -I$(srcdir)/../include @AM_CPPFLAGS@

noinst_LTLIBRARIES = libfactories.la

libfactories_la_SOURCES = \
	ArrayFactories.c \
	ArraySequenceFactories.c \
	StringVector.c \
	VectorFactories.c \
	VectorSequenceFactories.c

noinst_HEADERS = \
	CreateArraySequence_source.c \
	CreateArray_source.c \
	CreateVectorSequence_source.c \
	CreateVector_source.c \
	DestroyArraySequence_source.c \
	DestroyArray_source.c \
	DestroyVectorSequence_source.c \
	DestroyVector_source.c \
	ResizeArray_source.c \
	ResizeVector_source.c
