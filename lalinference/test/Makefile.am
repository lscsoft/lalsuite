AM_CFLAGS = $(LALSUPPORT_CFLAGS)
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalinference.la

if LALSUPPORT
if LALINSPIRAL
if LALXML

XMLPRG = LALInferenceXMLTest
endif

check_PROGRAMS = \
	InferenceTest $(XMLPRG)


TESTS = $(check_PROGRAMS) $(XML_PROGRAMS)
endif
endif

EXTRA_DIST =

CLEANFILES = \
	*.out \
	*.dat \
    *.xml

TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null
