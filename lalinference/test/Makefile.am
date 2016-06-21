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
if LALXML
test_programs += LALInferenceXMLTest
endif

# Add shell, Python, etc. test scripts to this variable
# Disable test_multiband.sh for now
#test_scripts = test_multiband.sh
if SWIG_BUILD_PYTHON
#test_scripts += test_plot.py
test_scripts += test_detframe.py
test_scripts += test_bayestar_sample_model_psd.sh
test_scripts += test_bayestar_realize_coincs.sh
test_scripts += test_bayestar_sim_to_tmpltbank.sh
# test_scripts += test_bayestar_lattice_tmpltbank.sh # FIXME: broken
if HAVE_CHEALPIX
test_scripts += test_bayestar.py
endif
endif

# Add any helper programs required by tests to this variable
test_helpers += LALInferenceMultiBandTest

# Note: This test file is included in the repository because lalapps cannot be a dependency of lalinference.
HL-INJECTIONS_1_TEST-1000000000-10.xml:
	lalapps_inspinj --m-distr fixMasses --fixed-mass1 1.4 --fixed-mass2 1.4 \
	--gps-start-time 1000000000 --gps-end-time 1000000010 --d-distr uniform \
	--min-distance 1e3 --max-distance 1e6 --l-distr random --i-distr uniform \
	--f-lower 30 --waveform TaylorF2threePointFivePN --disable-spin \
	--time-step 1 --user-tag TEST

test_bayestar_sample_model_psd.xml: test_bayestar_sample_model_psd.sh

test_bayestar_realize_coincs.sh: HL-INJECTIONS_1_TEST-1000000000-10.xml test_bayestar_sample_model_psd.xml

test_bayestar_realize_coincs.xml: test_bayestar_realize_coincs.sh

test_bayestar_sim_to_tmpltbank.sh: HL-INJECTIONS_1_TEST-1000000000-10.xml

MOSTLYCLEANFILES = \
	*.dat \
	*.out \
	*.xml \
	test_bayestar_sample_model_psd.xml \
	test_bayestar_realize_coincs.xml \
	test_bayestar_sim_to_tmpltbank.xml \
	test_bayestar_lattice_tmpltbank.xml \
	test_result_images/pp_plot_default.png \
	test_result_images/pp_plot_lines.png \
	test_result_images/pp_plot_steps.png \
	test_vot.xml \
	$(END_OF_LIST)

EXTRA_DIST += \
	LALInferenceTest.h \
	baseline_images/pp_plot_default.png \
	baseline_images/pp_plot_lines.png \
	baseline_images/pp_plot_steps.png \
	HL-INJECTIONS_1_TEST-1000000000-10.xml \
	$(END_OF_LIST)
