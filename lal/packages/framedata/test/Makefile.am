## Process this file with automake to produce Makefile.in
if FRAME
check_PROGRAMS   = FrameStreamTest Read40mData ReadGEOData ReadLHOData MakeFrames
TESTS            = FrameStreamTest
else
check_PROGRAMS   =
TESTS            = NoTest
endif

TESTS_ENVIRONMENT = LAL_FRAME_PATH=$(srcdir)

noinst_DATA = F-600000000-6-10.F F-600000060-6-10.F F-600000120-6-10.F

EXTRA_DIST = NoTest FrameDataTest.c SpecBufferTest.c DataBufferTest.c $(DATA)

LDADD = $(top_builddir)/packages/framedata/src/liblalframe.la $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

FrameStreamTest_SOURCES = FrameStreamTest.c
Read40mData_SOURCES = Read40mData.c
ReadGEOData_SOURCES = ReadGEOData.c
ReadLHOData_SOURCES = ReadLHOData.c
MakeFrames_SOURCES = MakeFrames.c

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep *.out *.[0-9][0-9][0-9]
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
