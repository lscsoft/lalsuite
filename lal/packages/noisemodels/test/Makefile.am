## Process this file with automake to produce Makefile.in
check_PROGRAMS = ChirpSpace RandomSignal 
ChirpSpace_SOURCES = ChirpSpace.c
ChirpSpace_LDADD = $(top_builddir)/lib/liblal.la
RandomSignal_SOURCES = RandomSignal.c
RandomSignal_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = ChirpSpace RandomSignal
EXTRA_DIST = RandomSignal.par demo 
CLEANFILES = *.out
