## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libfactories.la

libfactories_la_SOURCES = \
  VectorFactories.c VectorSequenceFactories.c \
  ArrayFactories.c ArraySequenceFactories.c
MAINTAINERCLEANFILES = $(libfactories_la_SOURCES)
EXTRA_DIST = \
  VectorFactories.m4 CreateVector.m4 ResizeVector.m4 DestroyVector.m4 \
  ArrayFactories.m4 CreateArray.m4 ResizeArray.m4 DestroyArray.m4 \
  VectorSequenceFactories.m4 CreateVectorSequence.m4 DestroyVectorSequence.m4 \
  ArraySequenceFactories.m4 CreateArraySequence.m4 DestroyArraySequence.m4

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep


VectorFactories.c: VectorFactories.m4 CreateVector.m4 DestroyVector.m4 ResizeVector.m4
	cd $(srcdir) && @M4@ VectorFactories.m4 > VectorFactories.c

VectorSequenceFactories.c: VectorSequenceFactories.m4 CreateVectorSequence.m4 DestroyVectorSequence.m4
	cd $(srcdir) && @M4@ VectorSequenceFactories.m4 > VectorSequenceFactories.c

ArrayFactories.c: ArrayFactories.m4 CreateArray.m4 DestroyArray.m4 ResizeArray.m4
	cd $(srcdir) && @M4@ ArrayFactories.m4 > ArrayFactories.c

ArraySequenceFactories.c: ArraySequenceFactories.m4 CreateArraySequence.m4 DestroyArraySequence.m4
	cd $(srcdir) && @M4@ ArraySequenceFactories.m4 > ArraySequenceFactories.c
@DVIMKFILE@
