## Process this file with automake to produce Makefile.in
check_PROGRAMS = LALHelloTest
LALHelloTest_SOURCES = LALHelloTest.c
LALHelloTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = LALHelloTest
TESTS_ENVIRONMENT = 1>/dev/null
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
