include $(top_srcdir)/gnuscripts/lalsuite_header_links.am

pkginclude_HEADERS = \
	AVFactories.h \
	SeqFactories.h \
	StringVector.h \
	$(END_OF_LIST)

noinst_LTLIBRARIES = libfactories.la

libfactories_la_SOURCES = \
	ArrayFactories.c \
	ArraySequenceFactories.c \
	StringVector.c \
	VectorFactories.c \
	VectorSequenceFactories.c \
	$(END_OF_LIST)

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
	ResizeVector_source.c \
	$(END_OF_LIST)
