## Process this file with automake to produce Makefile.in
check_PROGRAMS = InjectTest
EXTRA_DIST = StreamInput.h StreamVectorInput.c StreamVectorSequenceInput.c
InjectTest_SOURCES = InjectTest.c StreamVectorInput.c StreamVectorSequenceInput.c
InjectTest_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = InjectTest
