SUBDIRS = baseline_images

TEST_CPPFLAGS = -DTEST_DATA_DIR='"$(abs_srcdir)/"'
LDADD = $(top_builddir)/src/liblalinference.la

if LALXML

XMLPRG = \
	LALInferenceXMLTest

endif

LALInferenceTest_SOURCES = \
	LALInferenceTest.c \
	LALInferenceTest.h

LALInferencePriorTest_SOURCES = \
	LALInferencePriorTest.c
	LALInferenceTest.h

LALInferenceLikelihoodTest_SOURCES = \
	LALInferenceLikelihoodTest.c
	LALInferenceTest.h

LALInferenceProposalTest_SOURCES = \
	LALInferenceProposalTest.c
	LALInferenceTest.h

LALInferenceInjectionTest_SOURCES = \
	LALInferenceTestInjectionLikelihood.c
	LALInferenceTest.h

LALInferenceKDTest_SOURCES = \
	LALInferenceKDTest.c \
	LALInferenceTest.h

check_PROGRAMS = \
	LALInferenceTest \
	LALInferencePriorTest \
	LALInferenceLikelihoodTest \
	LALInferenceProposalTest \
	LALInferenceKDTest \
	LALInferenceInjectionTest \
	$(XMLPRG)

TESTS = \
	LALInferenceTest \
	LALInferencePriorTest \
	LALInferenceKDTest \
	$(XMLPRG)

EXTRA_DIST = test_plot.py

# FIXME: Run this as a regular unit test.
installcheck-local:
	PYTHONPATH=$(pyexecdir):$(PYTHONPATH) $(PYTHON) $(srcdir)/test_plot.py

CLEANFILES = \
	test_vot.xml \
	test_result_images/pp_plot_default.png \
	test_result_images/pp_plot_lines.png \
	test_result_images/pp_plot_steps.png \
	*.out \
	*.dat \
	*.xml
