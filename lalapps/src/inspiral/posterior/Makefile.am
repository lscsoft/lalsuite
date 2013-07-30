SUBDIRS = SPINspiral mpi

LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if LALINSPIRAL
if LALFRAME

if LALMETAIO
metaio_progs = \
	lalapps_inspnest
lalapps_inspnest_SOURCES = \
	lalnest.c \
	LALInspiralMCMC.c \
	LALInspiralMCMCUser.c \
	nest_calc.c \
	nest_calc.h \
	LALInspiralMCMC.h \
	LALInspiralMCMCUser.h
endif

if LALPULSAR
if LALINFERENCE
pulsar_progs = lalinference_nest ChiSquare_test #lalinference_mcmcmpi
ChiSquare_test_SOURCES = ChiSquare_test.c
lalinference_nest_SOURCES = LALInferenceNest.c

if HAVE_BAMBI
pulsar_progs += lalinference_bambi
lalinference_bambi_SOURCES = LALInferenceBAMBI.c
lalinference_bambi_LDADD = ${BAMBI_LIBS}
endif

endif

endif

frame_progs = \
	lalapps_followupMcmc
lalapps_followupMcmc_SOURCES = \
	followupMcmc.c

endif
endif


bin_PROGRAMS = $(metaio_progs) $(frame_progs) $(pulsar_progs) $(mpi_progs)

noinst_PROGRAMS = $(td_test_progs)

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_builddir)/src/lalapps -DPKG_DATA_DIR='"$(pkgdatadir)/"'

if HAVE_PYTHON
bin_SCRIPTS = \
	lalapps_combine_evidence \
	lalapps_coherence_test \
	lalapps_nest_pipe \
	lalapps_vo_nest2pos \
	lalinference_pipe \
	lalapps_nest2pos \
	lalapps_merge_nested_sampling_runs
pkgpython_DATA = \
	nest_utils.py \
	combine_evidence.py \
	nest2pos.py \
	lalinference_pipe_utils.py
endif

EXTRA_DIST = \
	LALInferenceMCMCSampler.c \
	LALInferenceMCMCSampler.h \
	combine_evidence.in \
	coherence_test.in \
	nest_pipe.in \
	nest_utils.py \
	combine_evidence.py \
	lalapps_vo_nest2pos.in \
	lalinference_pipe.in \
	lalinference_pipe_utils.py \
	nest2pos.in \
	nest2pos.py \
	merge_nested_sampling_runs.in

python_config.sed: Makefile
	@-rm -f python_config.sed
	@echo "s+@PYTHONLIBDIR@+@pkgpythondir@+g" > python_config.sed
	@echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed

lalapps_combine_evidence: $(srcdir)/combine_evidence.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/combine_evidence.in > lalapps_combine_evidence
	@chmod a+x lalapps_combine_evidence

lalapps_coherence_test: $(srcdir)/coherence_test.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/coherence_test.in > lalapps_coherence_test
	@chmod a+x lalapps_coherence_test

lalapps_nest_pipe: $(srcdir)/nest_pipe.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/nest_pipe.in >lalapps_nest_pipe
	@chmod a+x lalapps_nest_pipe

lalapps_vo_nest2pos: $(srcdir)/lalapps_vo_nest2pos.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/lalapps_vo_nest2pos.in >lalapps_vo_nest2pos
	@chmod a+x lalapps_vo_nest2pos

lalinference_pipe: $(srcdir)/lalinference_pipe.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/lalinference_pipe.in > lalinference_pipe
	@chmod a+x lalinference_pipe

lalapps_nest2pos: $(srcdir)/nest2pos.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/nest2pos.in > lalapps_nest2pos
	@chmod a+x lalapps_nest2pos

lalapps_merge_nested_sampling_runs: $(srcdir)/merge_nested_sampling_runs.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/merge_nested_sampling_runs.in >lalapps_merge_nested_sampling_runs
	@chmod a+x lalapps_merge_nested_sampling_runs


CLEANFILES = \
	$(bin_SCRIPTS) \
	python_config.sed

BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common

