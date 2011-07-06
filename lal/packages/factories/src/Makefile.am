noinst_LTLIBRARIES = libfactories.la

libfactories_la_SOURCES = \
	ArrayFactories.c \
	ArraySequenceFactories.c \
	StringVector.c \
	VectorFactories.c \
	VectorSequenceFactories.c

noinst_HEADERS = \
	CreateArray_source.c \
	DestroyArray_source.c \
	ResizeArray_source.c

MAINTAINERCLEANFILES = \
	ArraySequenceFactories.c \
	VectorFactories.c \
	VectorSequenceFactories.c

EXTRA_DIST = \
	ArraySequenceFactories.m4 \
	CreateArraySequence.m4 \
	DestroyArraySequence.m4 \
	VectorFactories.m4 \
	CreateVector.m4 \
	ResizeVector.m4 \
	DestroyVector.m4 \
	VectorSequenceFactories.m4 \
	CreateVectorSequence.m4 \
	DestroyVectorSequence.m4

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

VectorFactories.c: VectorFactories.m4 CreateVector.m4 DestroyVector.m4 ResizeVector.m4
	cd $(srcdir) && @M4@ VectorFactories.m4 > VectorFactories.c

VectorSequenceFactories.c: VectorSequenceFactories.m4 CreateVectorSequence.m4 DestroyVectorSequence.m4
	cd $(srcdir) && @M4@ VectorSequenceFactories.m4 > VectorSequenceFactories.c

ArraySequenceFactories.c: ArraySequenceFactories.m4 CreateArraySequence.m4 DestroyArraySequence.m4
	cd $(srcdir) && @M4@ ArraySequenceFactories.m4 > ArraySequenceFactories.c

@DVIMKFILE@
