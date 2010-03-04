AMCFLAGS = $(LALSUPPORT_CFLAGS)
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalinspiral.la

if LALSUPPORT
check_PROGRAMS = \
	BasicInjectTest \
	GeneratePPNInspiralTest \
	GeneratePPNAmpCorInspiralTest \
	RandomInspiralSignalTest
TESTS = $(check_PROGRAMS)
else
check_PROGRAMS =
TESTS =
endif

EXTRA_DIST = RandomInspiralSignalTest.par

CLEANFILES = *.out

TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null
