## Process this file with automake to produce Makefile.in
if MPI
check_PROGRAMS  = FindChirpTmpltTest FindChirpSPTest
TESTS           = FindChirpTmpltTest FindChirpSPTest
else
check_PROGRAMS  =
TESTS           = NoTest
endif

EXTRA_DIST      = NoTest FindChirpTmpltTest.c FindChirpSPTest.c

FindChirpTmpltTest_SOURCES = FindChirpTmpltTest.c
FindChirpTmpltTest_LDADD = $(top_builddir)/lib/liblal.la
FindChirpSPTest_SOURCES = FindChirpSPTest.c
FindChirpSPTest_LDADD = $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
CLEANFILES = FindChirpTmpltTest

