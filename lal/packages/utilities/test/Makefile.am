## Process this file with automake to produce Makefile.in
noinst_PROGRAMS = RandomTest FindRootTest IntegrateTest InterpolateTest SortTest
RandomTest_SOURCES = RandomTest.c
RandomTest_LDADD = $(top_srcdir)/lib/liblal.a
FindRootTest_SOURCES = FindRootTest.c
FindRootTest_LDADD = $(top_srcdir)/lib/liblal.a
IntegrateTest_SOURCES = IntegrateTest.c
IntegrateTest_LDADD = $(top_srcdir)/lib/liblal.a
InterpolateTest_SOURCES = InterpolateTest.c
InterpolateTest_LDADD = $(top_srcdir)/lib/liblal.a
SortTest_SOURCES = SortTest.c
SortTest_LDADD = $(top_srcdir)/lib/liblal.a
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
TESTS = $(PROGRAMS)
MOSTLYCLEANFILES = PrintVector.*
