## Process this file with automake to produce Makefile.in
check_PROGRAMS = ChirpSpace CoarseTest CoarseTest2 RandomSignal
ChirpSpace_SOURCES = ChirpSpace.c
ChirpSpace_LDADD = $(top_builddir)/lib/liblal.la
CoarseTest_SOURCES = CoarseTest.c
CoarseTest_LDADD = $(top_builddir)/lib/liblal.la
CoarseTest2_SOURCES = CoarseTest2.c
CoarseTest2_LDADD = $(top_builddir)/lib/liblal.la
RandomSignal_SOURCES = RandomSignal.c
RandomSignal_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = ChirpSpace CoarseTest CoarseTest2 RandomSignal
EXTRA_DIST = CoarseTest.par CoarseTest2.par demo demo2
CLEANFILES = *.out
