## Process this file with automake to produce Makefile.in
check_PROGRAMS = ChirpSpace CoarseTest CoarseTest2 MetricTest MetricTestBCV BCVTemplates PNTemplates BCVTemplatesFlatMesh InspiralSpinBankTest
ChirpSpace_SOURCES = ChirpSpace.c
ChirpSpace_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

CoarseTest_SOURCES = CoarseTest.c
CoarseTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

CoarseTest2_SOURCES = CoarseTest2.c
CoarseTest2_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

MetricTest_SOURCES = MetricTest.c
MetricTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

MetricTestBCV_SOURCES = MetricTestBCV.c
MetricTestBCV_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

BCVTemplates_SOURCES = BCVTemplates.c
BCVTemplates_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

PNTemplates_SOURCES = PNTemplates.c
PNTemplates_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

BCVTemplatesFlatMesh_SOURCES = BCVTemplatesFlatMesh.c
BCVTemplatesFlatMesh_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

InspiralSpinBankTest_SOURCES = InspiralSpinBankTest.c
if METAIO
InspiralSpinBankTest_LDADD = $(top_builddir)/packages/support/src/liblalmetaio.la $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
else
InspiralSpinBankTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
endif

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
#TESTS = ChirpSpace CoarseTest CoarseTest2 MetricTest MetricTestBCV PNTemplates BCVTemplates BCVTemplatesFlatMesh InspiralSpinBankTest
TESTS = ChirpSpace CoarseTest CoarseTest2 MetricTest MetricTestBCV PNTemplates BCVTemplates BCVTemplatesFlatMesh
EXTRA_DIST = CoarseTest.par demo 
MOSTLYCLEANFILES = .dvi-dep *.out
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
