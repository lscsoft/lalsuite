## Process this file with automake to produce Makefile.in
check_PROGRAMS = RandomTest MersenneRandomTest FindRootTest IntegrateTest InterpolateTest SortTest ODETest
RandomTest_SOURCES = RandomTest.c
RandomTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
MersenneRandomTest_SOURCES = MersenneRandomTest.c
MersenneRandomTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
FindRootTest_SOURCES = FindRootTest.c
FindRootTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
IntegrateTest_SOURCES = IntegrateTest.c
IntegrateTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
InterpolateTest_SOURCES = InterpolateTest.c
InterpolateTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
SortTest_SOURCES = SortTest.c
SortTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
ODETest_SOURCES = ODETest.c
ODETest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = RandomTest MersenneRandomTest FindRootTest IntegrateTest InterpolateTest SortTest ODETest
MOSTLYCLEANFILES = .dvi-dep PrintVector.* *.out
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
