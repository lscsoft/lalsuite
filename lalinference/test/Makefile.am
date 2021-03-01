EXTRA_DIST =
SUBDIRS =
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
test_programs += test_cubic_interp

# Add shell, Python, etc. test scripts to this variable
# Disable test_multiband.sh for now
# test_scripts = test_multiband.sh

# test lalinference in a higher level rather than unit tests

if SWIG_BUILD_PYTHON
SUBDIRS += python
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
