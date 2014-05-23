EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_test.am

# Add compiled test programs to this variable
test_programs += LALInferenceTest
test_programs += LALInferenceKDTest
test_programs += LALInferencePriorTest
#test_programs += LALInferenceInjectionTest
#test_programs += LALInferenceLikelihoodTest
#test_programs += LALInferenceProposalTest
if LALXML
test_programs += LALInferenceXMLTest
endif

# Add shell, Python, etc. test scripts to this variable
#test_scripts += test_plot.py
if HAVE_CHEALPIX
test_scripts += test_bayestar.py
endif

# Add any helper programs required by tests to this variable
test_helpers +=

MOSTLYCLEANFILES = \
	*.dat \
	*.out \
	*.xml \
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
	$(END_OF_LIST)
