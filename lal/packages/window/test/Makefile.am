## Process this file with automake to produce Makefile.in
check_PROGRAMS = WindowTest
WindowTest_SOURCES = WindowTest.c
WindowTest_LDADD = $(top_builddir)/lib/liblal.la -lm $(top_builddir)/packages/support/src/liblalsupport.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = WindowTest
MOSTLYCLEANFILES = .dvi-dep PrintVector.*
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
