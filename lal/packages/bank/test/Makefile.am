## Process this file with automake to produce Makefile.in
check_PROGRAMS = ChirpSpace CoarseTest 
ChirpSpace_SOURCES = ChirpSpace.c
ChirpSpace_LDADD = $(top_builddir)/lib/liblal.la
CoarseTest_SOURCES = CoarseTest.c
CoarseTest_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = ChirpSpace CoarseTest 
EXTRA_DIST = CoarseTest.par demo 
CLEANFILES = *.out
