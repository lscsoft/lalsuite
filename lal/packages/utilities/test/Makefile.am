## Process this file with automake to produce Makefile.in
check_PROGRAMS = RandomTest FindRootTest IntegrateTest InterpolateTest SortTest
RandomTest_SOURCES = RandomTest.c
RandomTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
FindRootTest_SOURCES = FindRootTest.c
FindRootTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
IntegrateTest_SOURCES = IntegrateTest.c
IntegrateTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
InterpolateTest_SOURCES = InterpolateTest.c
InterpolateTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
SortTest_SOURCES = SortTest.c
SortTest_LDADD = $(top_builddir)/lib/liblal.la $(top_builddir)/packages/support/src/liblalsupport.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = RandomTest FindRootTest IntegrateTest InterpolateTest SortTest
MOSTLYCLEANFILES = .dvi-dep PrintVector.*
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
