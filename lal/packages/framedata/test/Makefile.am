## Process this file with automake to produce Makefile.in
if FRAME
noinst_PROGRAMS = FrameDataTest SpecBufferTest DataBufferTest
TESTS           = $(PROGRAMS)
else
noinst_PROGRAMS =
TESTS           = NoTest
endif

EXTRA_DIST      = NoTest FrameDataTest.c SpecBufferTest.c DataBufferTest.c

FrameDataTest_SOURCES = FrameDataTest.c
SpecBufferTest_SOURCES = SpecBufferTest.c
DataBufferTest_SOURCES = DataBufferTest.c
FrameDataTest_LDADD = $(top_srcdir)/lib/liblal.a
SpecBufferTest_LDADD = $(top_srcdir)/lib/liblal.a
DataBufferTest_LDADD = $(top_srcdir)/lib/liblal.a
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = Segment.* Spectrum.* Response.*
