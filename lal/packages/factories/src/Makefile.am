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
	CreateVector_source.c \
	DestroyArraySequence_source.c \
	DestroyArray_source.c \
	DestroyVector_source.c \
	ResizeArray_source.c \
	ResizeVector_source.c

MAINTAINERCLEANFILES = \
	VectorSequenceFactories.c

EXTRA_DIST = \
	VectorSequenceFactories.m4 \
	CreateVectorSequence.m4 \
	DestroyVectorSequence.m4

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

VectorSequenceFactories.c: VectorSequenceFactories.m4 CreateVectorSequence.m4 DestroyVectorSequence.m4
	cd $(srcdir) && @M4@ VectorSequenceFactories.m4 > VectorSequenceFactories.c

@DVIMKFILE@
