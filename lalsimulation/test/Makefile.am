TEST_CPPFLAGS = -DTEST_DATA_DIR='"$(abs_srcdir)/"'
LDADD = $(top_builddir)/src/liblalsimulation.la

if HAVE_OPENMP
OPENMPPROGRAMS = TestOpenMP
endif

check_PROGRAMS = \
	$(OPENMPPROGRAMS) \
	GenerateSimulation \
	TestEOBNRv2 \
	GRFlagsTest \
	WaveformFlagsTest \
	TestSphHarmTS \
	TestWaveformFromCache \
	TestPrecessWaveform \
	TestPrecessWaveformIMRPhenomB \
	TestPrecessWaveformEOBNR \
	TestPhenomP \
	TestXLALSimAddInjection
TESTS = $(check_PROGRAMS)

# disable test programs
# TestSpinTaylorT4Dynamics - memory issues?

EXTRA_DIST =

MOSTLYCLEANFILES =
CLEANFILES = \
	*.dat \
	h_ref.txt \
	h_rot.txt \
	h_ref_EOBNR.txt \
	h_rot_EOBNR.txt \
	h_ref_PhenomB.txt \
	h_rot_PhenomB.txt
