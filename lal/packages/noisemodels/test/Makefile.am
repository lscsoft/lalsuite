## Process this file with automake to produce Makefile.in
check_PROGRAMS = NoisePSDTest
NoisePSDTest_SOURCES = NoisePSDTest.c
NoisePSDTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include -I/usr/local/include
TESTS = NoisePSDTest
TESTS_ENVIRONMENT=1>/dev/null 2>/dev/null
EXTRA_DIST = NoisePSDTest.eps
EXTRA_TEXSOURCES = NoisePSDTest.eps
MOSTLYCLEANFILES = .dvi-dep
CLEANFILES = *.out
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
