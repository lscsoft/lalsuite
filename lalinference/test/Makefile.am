AM_CFLAGS = $(LALSUPPORT_CFLAGS)
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalinspiral.la

if LALSUPPORT
if LALINSPIRAL
check_PROGRAMS = \
	InferenceTest
TESTS = $(check_PROGRAMS)
endif
endif

EXTRA_DIST =

CLEANFILES = \
	*.out \
	*.dat

TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null
