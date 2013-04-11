AM_CFLAGS = $(LALSUPPORT_CFLAGS)
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalsimulation.la

check_PROGRAMS = \
	GenerateSimulation \
	TestSpinTaylorT4Dynamics \
	TestEOBNRv2 \
	GRFlagsTest \
	WaveformFlagsTest \
	TestSphHarmTS \
	TestWaveformFromCache \
	TestPrecessWaveform
TESTS = $(check_PROGRAMS)

EXTRA_DIST =

MOSTLYCLEANFILES =
CLEANFILES = \
	*.dat \
	h_ref.txt \
	h_rot.txt
