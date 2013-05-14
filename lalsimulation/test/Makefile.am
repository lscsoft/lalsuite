AM_CFLAGS = $(LALSUPPORT_CFLAGS)
AM_CPPFLAGS = -I$(top_builddir)/include
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalsimulation.la

check_PROGRAMS = \
	GenerateSimulation \
	TestEOBNRv2 \
	GRFlagsTest \
	WaveformFlagsTest \
	TestSphHarmTS \
	TestWaveformFromCache \
	TestPrecessWaveformEOBNR
TESTS = $(check_PROGRAMS)

# disable test programs
# TestSpinTaylorT4Dynamics - memory issues?

EXTRA_DIST =

MOSTLYCLEANFILES =
CLEANFILES = \
	*.dat \
	h_ref.txt \
	h_rot.txt
