## Process this file with automake to produce Makefile.in
if FRAME
#check_PROGRAMS = FrameDataTest SpecBufferTest DataBufferTest
#TESTS          = FrameDataTest SpecBufferTest DataBufferTest
check_PROGRAMS   = FrameStreamTest
TESTS            = FrameStreamTest
else
check_PROGRAMS   =
TESTS            = NoTest
endif

TESTS_ENVIRONMENT = LAL_FRAME_PATH=$(srcdir)

noinst_DATA = P-682830700-2-6.F P-682830712-2-6.F P-682830724-2-6.F

EXTRA_DIST = NoTest FrameStreamTest.c FrameDataTest.c SpecBufferTest.c DataBufferTest.c $(DATA)

FrameStreamTest_SOURCES = FrameStreamTest.c
FrameStreamTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la

#FrameDataTest_SOURCES = FrameDataTest.c
#SpecBufferTest_SOURCES = SpecBufferTest.c
#DataBufferTest_SOURCES = DataBufferTest.c
#FrameDataTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
#SpecBufferTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
#DataBufferTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
#MOSTLYCLEANFILES = .dvi-dep Segment.* Spectrum.* Response.*
MOSTLYCLEANFILES = .dvi-dep *.out
#CLEANFILES = FrameDataTest SpecBufferTest DataBufferTest
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
