## Process this file with automake to produce Makefile.in
check_PROGRAMS = ChirpSpace FilterTest NoisePSDTest RandomSignal 
ChirpSpace_SOURCES = ChirpSpace.c
ChirpSpace_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
FilterTest_SOURCES = FilterTest.c
FilterTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
NoisePSDTest_SOURCES = NoisePSDTest.c
NoisePSDTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
RandomSignal_SOURCES = RandomSignal.c
RandomSignal_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = ChirpSpace FilterTest NoisePSDTest
EXTRA_DIST = RandomSignal.par demo
EXTRA_TEXSOURCES = FilterTest.c RandomSignal.c
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = *.out
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
