## Process this file with automake to produce Makefile.in
check_PROGRAMS = VectorOpsTest VectorIndexRangeTest
VectorOpsTest_SOURCES = VectorOpsTest.c
VectorOpsTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
VectorIndexRangeTest_SOURCES = VectorIndexRangeTest.c
VectorIndexRangeTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = VectorOpsTest VectorIndexRangeTest
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
