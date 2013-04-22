AM_CFLAGS = $(LALSUPPORT_CFLAGS)
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalsimulation.la

check_PROGRAMS = \
	GenerateSimulation \
	TestEOBNRv2 \
	GRFlagsTest \
	WaveformFlagsTest \
	TestSphHarmTS \
	TestWaveformFromCache \
	TestPrecessWaveform
TESTS = $(check_PROGRAMS)

# disable test programs
# TestSpinTaylorT4Dynamics - memory issues?

EXTRA_DIST =

MOSTLYCLEANFILES =
CLEANFILES = \
	*.dat \
	h_ref.txt \
	h_rot.txt
