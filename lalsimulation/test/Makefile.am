EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_test.am
AM_CPPFLAGS += -I$(top_srcdir)/lib

SUBDIRS = \
	$(END_OF_LIST)

# Add compiled test programs to this variable
test_programs += EOBNRv2Test
test_programs += GRFlagsTest
test_programs += LALSimulationTest
test_programs += PhenomPTest
test_programs += PhenomNSBHTest
test_programs += BHNSRemnantFitsTest
test_programs += NSBHPropertiesTest
test_programs += PNCoefficients
test_programs += PrecessWaveformEOBNRTest
test_programs += PrecessWaveformIMRPhenomBTest
test_programs += PrecessWaveformTest
test_programs += SphHarmTSTest
test_programs += WaveformFlagsTest
test_programs += WaveformFromCacheTest
test_programs += XLALSimAddInjectionTest
test_programs += InitialSpinRotationTest
test_programs += PrecessingHlmsTest
test_programs += SpinTaylorHlmsTest
test_programs += SEOBNRv4_ROM_NRTidalv2_NSBH_Test
#test_programs += TEOBResumROMTest
#test_programs += TestTaylorTFourier
#test_programs += SpinTaylorT4DynamicsTest
#if OPENMP
#test_programs += OpenMPTest
#endif

# Add shell, Python, etc. test scripts to this variable
# tests are currently broken on OS X
#test_scripts += AllApproximantsTest.sh
#if SWIG_BUILD_PYTHON
#test_scripts += LALSimWaveformReferenceTest.py
#endif

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

if HAVE_PYTHON
SUBDIRS += python
endif
