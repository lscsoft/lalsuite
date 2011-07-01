noinst_LTLIBRARIES = libfactories.la

libfactories_la_SOURCES = \
	ArrayFactories.c \
	ArraySequenceFactories.c \
	StringVector.c \
	VectorFactories.c \
	VectorSequenceFactories.c

MAINTAINERCLEANFILES = \
	ArrayFactories.c \
	ArraySequenceFactories.c \
	VectorFactories.c \
	VectorSequenceFactories.c

EXTRA_DIST = \
	ArrayFactories.m4 \
	CreateArray.m4 \
	ResizeArray.m4 \
	DestroyArray.m4 \
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

ArrayFactories.c: ArrayFactories.m4 CreateArray.m4 DestroyArray.m4 ResizeArray.m4
	cd $(srcdir) && @M4@ ArrayFactories.m4 > ArrayFactories.c

ArraySequenceFactories.c: ArraySequenceFactories.m4 CreateArraySequence.m4 DestroyArraySequence.m4
	cd $(srcdir) && @M4@ ArraySequenceFactories.m4 > ArraySequenceFactories.c

@DVIMKFILE@
