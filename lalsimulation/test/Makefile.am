AM_CFLAGS = $(LALSUPPORT_CFLAGS) @AM_CFLAGS@
AM_CPPFLAGS = -I$(top_builddir)/include -I$(top_srcdir)/src @AM_CPPFLAGS@
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalsimulation.la

check_PROGRAMS = \
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
