## Process this file with automake to produce Makefile.in
#check_PROGRAMS = FilterTest NoisePSDTest RandomSignal
check_PROGRAMS = NoisePSDTest
#FilterTest_SOURCES = FilterTest.c
NoisePSDTest_SOURCES = NoisePSDTest.c
#RandomSignal_SOURCES = RandomSignal.c
#FilterTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
NoisePSDTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
#RandomSignal_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
#TESTS = FilterTest NoisePSDTest RandomSignal
TESTS = NoisePSDTest
#EXTRA_DIST = RandomSignal.par
EXTRA_DIST = RandomSignal.par FilterTest.c RandomSignal.c
EXTRA_TEXSOURCES = FilterTest.c RandomSignal.c
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = *.out
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
