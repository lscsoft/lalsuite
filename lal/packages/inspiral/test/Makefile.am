## Process this file with automake to produce Makefile.in
check_PROGRAMS = LALInspiralTest ChirpSpace LALInspiralSpinningBHBinariesTest
LALInspiralTest_SOURCES = LALInspiralTest.c
ChirpSpace_SOURCES = ChirpSpace.c
LALInspiralSpinningBHBinariesTest_SOURCES =  LALInspiralSpinningBHBinariesTest.c
LALInspiralTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
ChirpSpace_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALInspiralSpinningBHBinariesTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = LALInspiralTest ChirpSpace LALInspiralSpinningBHBinariesTest
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = *.dat
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep

