## Process this file with automake to produce Makefile.in
noinst_PROGRAMS = OverlapTest DirichletTest CrossCorrTest
OverlapTest_SOURCES = OverlapTest.c
OverlapTest_LDADD = $(top_srcdir)/lib/liblal.a
DirichletTest_SOURCES = DirichletTest.c
DirichletTest_LDADD = $(top_srcdir)/lib/liblal.a
CrossCorrTest_SOURCES = CrossCorrTest.c
CrossCorrTest_LDADD = $(top_srcdir)/lib/liblal.a
INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
TESTS = $(PROGRAMS)
MOSTLYCLEANFILES = PrintVector.*
