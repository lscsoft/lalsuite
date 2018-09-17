EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_test.am

# Add compiled test programs to this variable
test_programs += LALInferenceTest
test_programs += LALInferencePriorTest
test_programs += LALInferenceGenerateROQTest
#test_programs += LALInferenceMultiBandTest
#test_programs += LALInferenceInjectionTest
#test_programs += LALInferenceLikelihoodTest
#test_programs += LALInferenceProposalTest
test_programs += LALInferenceHDF5Test

# Add shell, Python, etc. test scripts to this variable
# Disable test_multiband.sh for now
#test_scripts = test_multiband.sh
if SWIG_BUILD_PYTHON
test_scripts += test_detframe.py
endif

# Add any helper programs required by tests to this variable
test_helpers += LALInferenceMultiBandTest

MOSTLYCLEANFILES = \
	*.dat \
	*.out \
	test.hdf5 \
	$(END_OF_LIST)

EXTRA_DIST += \
	LALInferenceTest.h \
	$(END_OF_LIST)
