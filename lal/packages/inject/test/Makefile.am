## Process this file with automake to produce Makefile.in
check_PROGRAMS = BasicInjectTest InjectTest GeneratePPNInspiralTest SimulateTaylorCWTest SimulateSBTest
InjectTest_SOURCES = InjectTest.c
InjectTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
BasicInjectTest_SOURCES = BasicInjectTest.c
BasicInjectTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
GeneratePPNInspiralTest_SOURCES = GeneratePPNInspiralTest.c
GeneratePPNInspiralTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
SimulateTaylorCWTest_SOURCES = SimulateTaylorCWTest.c
SimulateTaylorCWTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
SimulateSBTest_SOURCES = SimulateSBTest.c CheckStatus.c
SimulateSBTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = BasicInjectTest GeneratePPNInspiralTest SimulateTaylorCWTest SimulateSBTest
EXTRA_DIST = CheckStatus.c CheckStatus.h
MOSTLYCLEANFILES = .dvi-dep WHITENED-SB*
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
