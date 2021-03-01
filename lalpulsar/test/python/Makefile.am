include $(top_srcdir)/gnuscripts/lalsuite_test_python.am

# Add compiled test programs to this variable
test_programs +=

# Add shell, Python, etc. test scripts to this variable
test_scripts += \
	test_computePSD.py \
	test_heterodyned_model.py \
	test_ReadSFDB.py \
	test_simulateCW.py \
	$(END_OF_LIST)

# Add any helper programs required by tests to this variable
test_helpers +=

MOSTLYCLEANFILES = \
	V-1_V1_1800SFT_simCW_simulateCWTest-*.sft \
	$(END_OF_LIST)

# Add references files required by tests to this variable
# to make them part of the distribution tarball
EXTRA_DIST += V1_mfdv5_20191114_043831.SFDB09
