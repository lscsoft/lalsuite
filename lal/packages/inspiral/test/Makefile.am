## Process this file with automake to produce Makefile.in
check_PROGRAMS = LALInspiralTest ChirpSpace
LALInspiralTest_SOURCES = LALInspiralTest.c
ChirpSpace_SOURCES = ChirpSpace.c
LALInspiralTest_LDADD = $(top_builddir)/lib/liblal.la
ChirpSpace_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = LALInspiralTest ChirpSpace
CLEANFILES = *.dat
