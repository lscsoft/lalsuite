## Process this file with automake to produce Makefile.in
#
# $Id$
#
check_PROGRAMS  = FindChirpBCVSpinTest FindChirpTDTest 
#CoherentInspiralInputTest
TESTS           = NoTest.sh

EXTRA_DIST      = NoTest.sh FindChirpTmpltTest.c FindChirpSPTest.c FindChirpBCVSpinTest.c 

EXTRA_TEXSOURCES = NoTest.sh

#FindChirpTmpltTest_SOURCES = FindChirpTmpltTest.c
#FindChirpTmpltTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
#CoherentInspiralInputTest_SOURCES = CoherentInspiralInputTest.c
FindChirpTDTest_SOURCES = FindChirpTDTest.c FindChirpTDTest.h
FindChirpTDTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
FindChirpBCVSpinTest_SOURCES = FindChirpBCVSpinTest.c
LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = FindChirpTmpltTest
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
