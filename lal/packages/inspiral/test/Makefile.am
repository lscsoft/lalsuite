## Process this file with automake to produce Makefile.in
check_PROGRAMS = GenerateInspiralWaveform ChirpSpace LALInspiralTest LALInspiralSpinningBHBinariesTest LALSTPNWaveformTest 

GenerateInspiralWaveform_SOURCES = GenerateInspiralWaveform.c
LALSTPNWaveformTest_SOURCES = LALSTPNWaveformTest.c
LALSTPNWaveformTest_LDADD= $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
ChirpSpace_SOURCES = ChirpSpace.c
LALInspiralSpinningBHBinariesTest_SOURCES =  LALInspiralSpinningBHBinariesTest.c 

LALInspiralTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
ChirpSpace_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
GenerateInspiralWaveform_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALInspiralSpinningBHBinariesTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
#TESTS = LALInspiralTest ChirpSpace LALInspiralTestOne LALInspiralSpinningBHBinariesTest LALSTPNWaveformTest 
TESTS = GenerateInspiralWaveform ChirpSpace LALInspiralTest LALInspiralSpinningBHBinariesTest
TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = *.dat
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep

@DVIMKFILE@
