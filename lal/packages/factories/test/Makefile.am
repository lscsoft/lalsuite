## Process this file with automake to produce Makefile.in
check_PROGRAMS = \
  VectorFactoriesTest \
  VectorSequenceFactoriesTest \
  ArrayFactoriesTest \
  ArraySequenceFactoriesTest

VectorFactoriesTest_SOURCES = VectorFactoriesTest.c
VectorSequenceFactoriesTest_SOURCES = VectorSequenceFactoriesTest.c

VectorFactoriesTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
VectorSequenceFactoriesTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

ArrayFactoriesTest_SOURCES = ArrayFactoriesTest.c
ArraySequenceFactoriesTest_SOURCES = ArraySequenceFactoriesTest.c

ArrayFactoriesTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
ArraySequenceFactoriesTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = VectorFactoriesTest \
        VectorSequenceFactoriesTest \
        ArrayFactoriesTest \
        ArraySequenceFactoriesTest
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep


EXTRA_DIST = \
        VectorFactoriesTest.m4 VectorFactoriesTestFunction.m4 \
        VectorSequenceFactoriesTest.m4 VectorSequenceFactoriesTestFunction.m4 \
        ArrayFactoriesTest.m4 ArrayFactoriesTestFunction.m4 \
        ArraySequenceFactoriesTest.m4 ArraySequenceFactoriesTestFunction.m4
MAINTAINERCLEANFILES = $(SOURCES)

VectorFactoriesTest.c: VectorFactoriesTest.m4 VectorFactoriesTestFunction.m4
	cd $(srcdir) && @M4@ VectorFactoriesTest.m4 > VectorFactoriesTest.c

VectorSequenceFactoriesTest.c: VectorSequenceFactoriesTest.m4 VectorSequenceFactoriesTestFunction.m4
	cd $(srcdir) && @M4@ VectorSequenceFactoriesTest.m4 > VectorSequenceFactoriesTest.c

ArrayFactoriesTest.c: ArrayFactoriesTest.m4 ArrayFactoriesTestFunction.m4
	cd $(srcdir) && @M4@ ArrayFactoriesTest.m4 > ArrayFactoriesTest.c

ArraySequenceFactoriesTest.c: ArraySequenceFactoriesTest.m4 ArraySequenceFactoriesTestFunction.m4
	cd $(srcdir) && @M4@ ArraySequenceFactoriesTest.m4 > ArraySequenceFactoriesTest.c
@DVIMKFILE@
