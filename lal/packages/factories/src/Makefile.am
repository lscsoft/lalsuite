## Process this file with automake to produce Makefile.in
noinst_LIBRARIES = liblal.a

liblal_a_SOURCES = \
  VectorFactories.c VectorSequenceFactories.c \
  ArrayFactories.c
MAINTAINERCLEANFILES = $(liblal_a_SOURCES)
EXTRA_DIST = \
  VectorFactories.m4 CreateVector.m4 DestroyVector.m4 \
  ArrayFactories.m4 CreateArray.m4 DestroyArray.m4 \
  VectorSequenceFactories.m4 CreateVectorSequence.m4 DestroyVectorSequence.m4

LIBLAL = $(top_srcdir)/lib/liblal.a
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include


.PHONY: liblal.a
liblal.a: $(liblal_a_OBJECTS) $(liblal_a_DEPENDENCIES)
	$(AR) cru $(LIBLAL) $(liblal_a_OBJECTS) $(liblal_a_LIBADD)
	$(RANLIB) $(LIBLAL)

VectorFactories.c: VectorFactories.m4 CreateVector.m4 DestroyVector.m4
	m4 VectorFactories.m4 > VectorFactories.c

VectorSequenceFactories.c: VectorSequenceFactories.m4 CreateVectorSequence.m4 DestroyVectorSequence.m4
	m4 VectorSequenceFactories.m4 > VectorSequenceFactories.c

ArrayFactories.c: ArrayFactories.m4 CreateArray.m4 DestroyArray.m4
	m4 ArrayFactories.m4 > ArrayFactories.c
