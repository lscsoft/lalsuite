## Process this file with automake to produce Makefile.in
check_PROGRAMS = ChirpSpace CoarseTest CoarseTest2
ChirpSpace_SOURCES = ChirpSpace.c
ChirpSpace_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
CoarseTest_SOURCES = CoarseTest.c
CoarseTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
CoarseTest2_SOURCES = CoarseTest2.c
CoarseTest2_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = ChirpSpace CoarseTest CoarseTest2
#EXTRA_DIST = CoarseTest2.par CoarseTest.par demo 
EXTRA_DIST = CoarseTest.par demo 
MOSTLYCLEANFILES = .dvi-dep *.out
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
