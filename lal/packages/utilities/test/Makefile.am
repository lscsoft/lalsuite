## Process this file with automake to produce Makefile.in
check_PROGRAMS = RandomTest FindRootTest IntegrateTest InterpolateTest SortTest
RandomTest_SOURCES = RandomTest.c
RandomTest_LDADD = $(top_srcdir)/lib/liblal.la
FindRootTest_SOURCES = FindRootTest.c
FindRootTest_LDADD = $(top_srcdir)/lib/liblal.la
IntegrateTest_SOURCES = IntegrateTest.c
IntegrateTest_LDADD = $(top_srcdir)/lib/liblal.la
InterpolateTest_SOURCES = InterpolateTest.c
InterpolateTest_LDADD = $(top_srcdir)/lib/liblal.la
SortTest_SOURCES = SortTest.c
SortTest_LDADD = $(top_srcdir)/lib/liblal.la
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
TESTS = RandomTest FindRootTest IntegrateTest InterpolateTest SortTest
MOSTLYCLEANFILES = PrintVector.*
