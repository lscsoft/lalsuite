TEST_CPPFLAGS = -DTEST_DATA_DIR='"$(abs_srcdir)/"'
LDADD = $(top_builddir)/src/liblalinference.la

if LALXML

XMLPRG = \
	LALInferenceXMLTest

endif

LALInferenceTest_SOURCES = \
	LALInferenceTest.c \
	LALInferenceTest.h

LALInferencePriorTest_SOURCES = \
	LALInferencePriorTest.c
	LALInferenceTest.h

LALInferenceLikelihoodTest_SOURCES = \
	LALInferenceLikelihoodTest.c
	LALInferenceTest.h

LALInferenceProposalTest_SOURCES = \
	LALInferenceProposalTest.c
	LALInferenceTest.h

LALInferenceInjectionTest_SOURCES = \
	LALInferenceTestInjectionLikelihood.c
	LALInferenceTest.h

LALInferenceKDTest_SOURCES = \
	LALInferenceKDTest.c \
	LALInferenceTest.h

check_PROGRAMS = \
	LALInferenceTest \
	LALInferencePriorTest \
	LALInferenceLikelihoodTest \
	LALInferenceProposalTest \
	LALInferenceKDTest \
	LALInferenceInjectionTest \
	$(XMLPRG)

TESTS = \
	LALInferenceTest \
	LALInferenceKDTest \
	$(XMLPRG)

#	LALInferencePriorTest

EXTRA_DIST =

CLEANFILES = \
	test_vot.xml \
	*.out \
	*.dat \
	*.xml
