## Process this file with automake to produce Makefile.in
check_PROGRAMS = InjectTest GeneratePPNInspiralTest
EXTRA_DIST = StreamInput.h StreamVectorInput.c StreamVectorSequenceInput.c
InjectTest_SOURCES = InjectTest.c StreamVectorInput.c StreamVectorSequenceInput.c
InjectTest_LDADD = $(top_builddir)/lib/liblal.la
GeneratePPNInspiralTest_SOURCES = GeneratePPNInspiralTest.c
GeneratePPNInspiralTest_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = InjectTest GeneratePPNInspiralTest
