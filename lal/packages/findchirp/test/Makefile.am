## Process this file with automake to produce Makefile.in
#
# $Id$
#
if MPI
check_PROGRAMS  = FindChirpBCVTest FindChirpBCVSpinTest
TESTS           = NoTest.sh
else
check_PROGRAMS  = FindChirpBCVTest FindChirpBCVSpinTest
TESTS           = NoTest.sh
endif

EXTRA_DIST      = NoTest.sh FindChirpTmpltTest.c FindChirpSPTest.c TwoInterfFindChirpSPTest.c FindChirpBCVTest.c FindChirpBCVSpinTest.c

EXTRA_TEXSOURCES = NoTest.sh

#FindChirpTmpltTest_SOURCES = FindChirpTmpltTest.c
#FindChirpTmpltTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
#FindChirpSPTest_SOURCES = FindChirpSPTest.c
#FindChirpSPTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
FindChirpBCVTest_SOURCES = FindChirpBCVTest.c
FindChirpBCVSpinTest_SOURCES = FindChirpBCVSpinTest.c
FindChirpBCVTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
FindChirpBCVSpinTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = FindChirpTmpltTest
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
