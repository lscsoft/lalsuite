## Process this file with automake to produce Makefile.in
#
# $Id$
#
if MPI
check_PROGRAMS  = FindChirpTmpltTest FindChirpSPTest
TESTS           = FindChirpTmpltTest FindChirpSPTest
else
check_PROGRAMS  =
TESTS           = NoTest
endif

EXTRA_DIST      = NoTest FindChirpTmpltTest.c FindChirpSPTest.c

FindChirpTmpltTest_SOURCES = FindChirpTmpltTest.c
FindChirpTmpltTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
FindChirpSPTest_SOURCES = FindChirpSPTest.c
FindChirpSPTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = FindChirpTmpltTest
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
