EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_test.am
AM_CPPFLAGS = -I$(top_srcdir)/src @AM_CPPFLAGS@ 

# Add compiled test programs to this variable
test_programs += EOBNRv2Test
test_programs += GRFlagsTest
test_programs += PhenomPTest
test_programs += PNCoefficients
test_programs += PrecessWaveformEOBNRTest
test_programs += PrecessWaveformIMRPhenomBTest
test_programs += PrecessWaveformTest
test_programs += SphHarmTSTest
test_programs += WaveformFlagsTest
test_programs += WaveformFromCacheTest
test_programs += XLALSimAddInjectionTest
#test_programs += SpinTaylorT4DynamicsTest
if HAVE_OPENMP
test_programs += OpenMPTest
endif

# Add shell, Python, etc. test scripts to this variable
test_scripts += AllApproximantsTest.sh
if SWIG_BUILD_PYTHON
test_scripts += LALSimWaveformReferenceTest.py
endif

# Add any helper programs required by tests to this variable
test_helpers += GenerateSimulation

MOSTLYCLEANFILES = \
	*.dat \
	h_ref.txt \
	h_ref_EOBNR.txt \
	h_ref_PhenomB.txt \
	h_rot.txt \
	h_rot_EOBNR.txt \
	h_rot_PhenomB.txt \
	$(END_OF_LIST)

EXTRA_DIST += \
	reviewed_waveforms.asc
	$(END_OF_LIST)
