AM_CFLAGS = $(LALSUPPORT_CFLAGS)
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalinspiral.la

if LALSUPPORT
check_PROGRAMS = \
	BasicInjectTest \
	GeneratePPNInspiralTest \
	GeneratePPNAmpCorInspiralTest \
	InjectionInterfaceTest \
	RandomInspiralSignalTest
TESTS = $(check_PROGRAMS)
else
check_PROGRAMS =
TESTS =
endif

EXTRA_DIST = \
	NRWaveIOTest.cfg \
	NRWaveIOTest.data \
	RandomInspiralSignalTest.par \
	injection.xml

CLEANFILES = \
	*.out \
	*.dat

TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null
