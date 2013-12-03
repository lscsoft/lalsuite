AM_CFLAGS = $(LALSUPPORT_CFLAGS)
AM_CPPFLAGS = -I$(top_builddir)/include
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalinference.la

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
	LALInferencePriorTest \
	LALInferenceKDTest \
	$(XMLPRG)

EXTRA_DIST =

CLEANFILES = \
	test_vot.xml \
	*.out \
	*.dat \
	*.xml
