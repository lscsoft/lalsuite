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
#test_scripts += test_plot.py
test_scripts += test_detframe.py
if HAVE_CHEALPIX
test_scripts += doctests.py

# installcheck-local:
# 	../python/bayestar_sample_model_psd --H1=aLIGOZeroDetHighPower --L1=aLIGOZeroDetHighPower -o test_bayestar_sample_model_psd.xml && \
# 	../python/bayestar_realize_coincs --detector H1 L1 --reference-psd test_bayestar_sample_model_psd.xml $(srcdir)/HL-INJECTIONS_1_TEST-1000000000-10.xml -o test_bayestar_realize_coincs.xml && \
# 	../python/bayestar_sim_to_tmpltbank $(srcdir)/HL-INJECTIONS_1_TEST-1000000000-10.xml -o test_bayestar_sim_to_tmpltbank.xml && \
# 	../python/bayestar_prune_neighborhood_tmpltbank test_bayestar_sim_to_tmpltbank.xml --mass1 1.4 --mass2 1.4 --snr 8 -o test_bayestar_prune_neighborhood_tmpltbank.xml
# bayestar_lattice_tmpltbank -o test_bayestar_lattice_tmpltbank.xml --min-mass 2.7 --max-mass 2.9
# bayestar_littlehope $srcdir/HL-INJECTIONS_1_TEST-1000000000-10.xml \
#     --template-bank test_bayestar_sim_to_tmpltbank.xml \
#     --reference-psd test_bayestar_sample_model_psd.xml \
#     --detector H1 --detector L1 --waveform TaylorF2threePointFivePN \
#     -o test_bayestar_littlehope.xml
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

MOSTLYCLEANFILES = \
	*.dat \
	*.out \
	test.hdf5 \
	test_bayestar_sample_model_psd.xml \
	test_bayestar_realize_coincs.xml \
	test_bayestar_sim_to_tmpltbank.xml \
	test_bayestar_prune_neighborhood_tmpltbank.xml \
	test_bayestar_littlehope.xml \
	test_bayestar_lattice_tmpltbank.xml \
	test_result_images/pp_plot_default.png \
	test_result_images/pp_plot_lines.png \
	test_result_images/pp_plot_steps.png \
	$(END_OF_LIST)

EXTRA_DIST += \
	LALInferenceTest.h \
	baseline_images/pp_plot_default.png \
	baseline_images/pp_plot_lines.png \
	baseline_images/pp_plot_steps.png \
	HL-INJECTIONS_1_TEST-1000000000-10.xml \
	$(END_OF_LIST)
