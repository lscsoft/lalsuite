## Process this file with automake to produce Makefile.in
if FRAME
#check_PROGRAMS = FrameDataTest SpecBufferTest DataBufferTest
#TESTS          = FrameDataTest SpecBufferTest DataBufferTest
check_PROGRAMS   = FrameStreamTest Read40mData ReadGEOData ReadLHOData MakeFrames
TESTS            = FrameStreamTest
else
check_PROGRAMS   =
TESTS            = NoTest
endif

TESTS_ENVIRONMENT = LAL_FRAME_PATH=$(srcdir)

noinst_DATA = F-600000000-6-10.F F-600000060-6-10.F F-600000120-6-10.F

EXTRA_DIST = NoTest FrameDataTest.c SpecBufferTest.c DataBufferTest.c $(DATA)

FrameStreamTest_SOURCES = FrameStreamTest.c
FrameStreamTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

Read40mData_SOURCES = Read40mData.c
Read40mData_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

ReadGEOData_SOURCES = ReadGEOData.c
ReadGEOData_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

ReadLHOData_SOURCES = ReadLHOData.c
ReadLHOData_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

MakeFrames_SOURCES = MakeFrames.c
MakeFrames_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

#FrameDataTest_SOURCES = FrameDataTest.c
#SpecBufferTest_SOURCES = SpecBufferTest.c
#DataBufferTest_SOURCES = DataBufferTest.c
#FrameDataTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
#SpecBufferTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
#DataBufferTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
#MOSTLYCLEANFILES = .dvi-dep Segment.* Spectrum.* Response.*
MOSTLYCLEANFILES = .dvi-dep *.out *.[0-9][0-9][0-9]
#CLEANFILES = FrameDataTest SpecBufferTest DataBufferTest
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
