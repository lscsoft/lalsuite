## Process this file with automake to produce Makefile.in
check_PROGRAMS = FilterTest NoisePSDTest RandomInspiralSignalTest 
FilterTest_SOURCES = FilterTest.c 
FilterTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
NoisePSDTest_SOURCES = NoisePSDTest.c
NoisePSDTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
RandomInspiralSignalTest_SOURCES = RandomInspiralSignalTest.c
RandomInspiralSignalTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include -I/usr/local/include
TESTS = FilterTest NoisePSDTest RandomInspiralSignalTest
#EXTRA_DIST = NoisePSDTest.par RandomInspiralSignalTest.par demo
EXTRA_DIST = RandomInspiralSignalTest.par demo NoisePSDTest.eps
EXTRA_TEXSOURCES = NoisePSDTest.eps
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = *.out
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
