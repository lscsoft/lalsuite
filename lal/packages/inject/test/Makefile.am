## Process this file with automake to produce Makefile.in
check_PROGRAMS = InjectTest BasicInjectTest GeneratePPNInspiralTest
InjectTest_SOURCES = InjectTest.c
InjectTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
BasicInjectTest_SOURCES = BasicInjectTest.c
BasicInjectTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
GeneratePPNInspiralTest_SOURCES = GeneratePPNInspiralTest.c
GeneratePPNInspiralTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = InjectTest BasicInjectTest GeneratePPNInspiralTest
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
