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

check_PROGRAMS = \
	LALInferenceTest \
	$(XMLPRG)

#	LALInferenceProposalTest \
#	LALInferenceKDTest

TESTS = \
	$(check_PROGRAMS) \
	$(XML_PROGRAMS)

endif
endif

EXTRA_DIST =

CLEANFILES = \
	test_vot.xml \
	*.out \
	*.dat \
	*.xml

TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null
