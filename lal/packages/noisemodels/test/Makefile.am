## Process this file with automake to produce Makefile.in
check_PROGRAMS = NoisePSDTest RandomSignal
NoisePSDTest_SOURCES = NoisePSDTest.c
RandomSignal_SOURCES = RandomSignal.c
NoisePSDTest_LDADD = $(top_builddir)/lib/liblal.la
RandomSignal_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = NoisePSDTest RandomSignal
EXTRA_DIST = RandomSignal.par
CLEANFILES = *.out
