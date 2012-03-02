AM_CFLAGS = $(LALSUPPORT_CFLAGS)
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalinference.la

if LALSUPPORT
if LALINSPIRAL
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

#LALInferenceProposalTest_SOURCES = \
#	LALInferenceProposalTest.c
#	LALInferenceTest.h

check_PROGRAMS = \
	LALInferenceTest \
	LALInferenceLikelihoodTest\
	$(XMLPRG) \
	LALInferencePriorTest
#	LALInferenceProposalTest \
#	LALInferenceKDTest

TESTS = \
	$(check_PROGRAMS)

endif
endif

EXTRA_DIST =

CLEANFILES = \
	test_vot.xml \
	*.out \
	*.dat \
	*.xml

TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null
