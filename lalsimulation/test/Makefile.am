AMCFLAGS = $(LALSUPPORT_CFLAGS)
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalsimulation.la

if LALSUPPORT
check_PROGRAMS = \
	GenerateInspiralWaveform
TESTS = $(check_PROGRAMS)
endif

EXTRA_DIST =

MOSTLYCLEANFILES =
CLEANFILES = \
	*.dat
