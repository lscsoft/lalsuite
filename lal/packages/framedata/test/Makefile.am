## Process this file with automake to produce Makefile.in
if FRAME
check_PROGRAMS = FrameDataTest SpecBufferTest DataBufferTest
TESTS          = FrameDataTest SpecBufferTest DataBufferTest
else
check_PROGRAMS =
TESTS           = NoTest
endif

EXTRA_DIST      = NoTest FrameDataTest.c SpecBufferTest.c DataBufferTest.c

FrameDataTest_SOURCES = FrameDataTest.c
SpecBufferTest_SOURCES = SpecBufferTest.c
DataBufferTest_SOURCES = DataBufferTest.c
FrameDataTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
SpecBufferTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
DataBufferTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep Segment.* Spectrum.* Response.*
CLEANFILES = FrameDataTest SpecBufferTest DataBufferTest
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
