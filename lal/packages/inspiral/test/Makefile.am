## Process this file with automake to produce Makefile.in
#if FRAME
if METAIO
#check_PROGRAMS = GenerateInspiralWaveform LALHybridTest
#else
#check_PROGRAMS = GenerateInspiralWaveform
#endif
#else
check_PROGRAMS = GenerateInspiralWaveform
endif

GenerateInspiralWaveform_SOURCES = GenerateInspiralWaveform.c
#LALHybridTest_SOURCES = LALHybridTest.c
GenerateInspiralWaveform_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
#LALHybridTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = GenerateInspiralWaveform
EXTRA_DIST = LALHybridTest.c
TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = *.dat
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep

@DVIMKFILE@
