## Process this file with automake to produce Makefile.in
check_PROGRAMS = BlockTest
BlockTest_SOURCES = BlockTest.c
BlockTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = BlockTest
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
