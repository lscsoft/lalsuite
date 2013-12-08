## Process this file with automake to produce Makefile.in

SUBDIRS = bayestar posterior

INSPIRALSRC = inspiral.c inspiralutils.c inspiral.h
TMPLTBANKSRC = tmpltbank.c inspiralutils.c inspiral.h
COHERESRC = coherent_inspiral.c
COHIRESRC = cohire.c
COHEREBANK = coherentbank.c
INCASRC = inca.c
THINCASRC = thinca.c
SIRESRC = sire.c
COIRESRC = coire.c
SIVASRC = siva.c
INSPINJSRC = inspinj.c inspiralutils.c inspiral.h
INJCUTSRC = injcut.c
MINJSRC = minj.c
BBHINJSRC = bbhinj.c
BLINDINJSRC = blindinj.c
SPLITBANKSRC = splitbank.c
TRIGBANKSRC  = trigbank.c
CHISQSRC = populate_chisq.c
TRIGSCANSRC = trigscan.c
INSPFRINJSRC = inspfrinj.c
GWF2XMLSRC = gwf2xml.c
SPININJSRC = spininj.c
COINEXTSRC = coinext.c
RANDOMBANKSRC = randombank.c
NINJASRC = ninja.c
FRNINJASRC = fr_ninja.c
MDCNINJASRC = mdc_ninja.c inspiralutils.c inspiral.h
CALCEXPSNRSRC = calcexpsnr.c
SNEDSRC = sned.c
CORSESRC = corse.c
SKYMAPSRC = skymap.c
COINJSRC = coinj.c
STOCHBANKSRC = stochasticbank.c

if LALFRAME
if LALMETAIO
if LALINSPIRAL
LALFRAME_LALMETAIO_PROGS = \
	lalapps_inspiral \
	lalapps_tmpltbank \
	lalapps_inca \
	lalapps_sire \
	lalapps_coire \
	lalapps_siva \
	lalapps_inspinj \
	lalapps_injcut \
	lalapps_splitbank \
	lalapps_minj \
	lalapps_bbhinj \
	lalapps_blindinj \
	lalapps_inspfrinj \
	lalapps_gwf2xml \
	lalapps_thinca \
	lalapps_trigbank \
	lalapps_populate_chisq \
	lalapps_trigscan \
	lalapps_spininj \
	lalapps_cohire \
	lalapps_coherentbank \
	lalapps_coinext \
	lalapps_randombank \
	lalapps_ninja \
	lalapps_fr_ninja \
	lalapps_mdc_ninja \
	lalapps_calcexpsnr \
	lalapps_sned \
	lalapps_corse \
	lalapps_skymap \
	lalapps_coinj \
	lalapps_cbc_stochasticbank
lalapps_inspiral_SOURCES = $(INSPIRALSRC)
lalapps_cohire_SOURCES = $(COHIRESRC)
lalapps_coherentbank_SOURCES = $(COHEREBANK)
lalapps_tmpltbank_SOURCES = $(TMPLTBANKSRC)
lalapps_inca_SOURCES = $(INCASRC)
lalapps_thinca_SOURCES = $(THINCASRC)
lalapps_sire_SOURCES = $(SIRESRC)
lalapps_coire_SOURCES = $(COIRESRC)
lalapps_siva_SOURCES = $(SIVASRC)
lalapps_inspinj_SOURCES = $(INSPINJSRC)
lalapps_injcut_SOURCES = $(INJCUTSRC)
lalapps_minj_SOURCES = $(MINJSRC)
lalapps_bbhinj_SOURCES = $(BBHINJSRC)
lalapps_blindinj_SOURCES = $(BLINDINJSRC)
lalapps_splitbank_SOURCES = $(SPLITBANKSRC)
lalapps_trigbank_SOURCES = $(TRIGBANKSRC)
lalapps_populate_chisq_SOURCES = $(CHISQSRC)
lalapps_trigscan_SOURCES = $(TRIGSCANSRC)
lalapps_spininj_SOURCES = $(SPININJSRC)
lalapps_inspfrinj_SOURCES = $(INSPFRINJSRC)
lalapps_gwf2xml_SOURCES = $(GWF2XMLSRC)
lalapps_coinext_SOURCES = $(COINEXTSRC)
lalapps_randombank_SOURCES = $(RANDOMBANKSRC)
lalapps_ninja_SOURCES = $(NINJASRC)
lalapps_fr_ninja_SOURCES = $(FRNINJASRC)
lalapps_mdc_ninja_SOURCES = $(MDCNINJASRC)
lalapps_calcexpsnr_SOURCES = $(CALCEXPSNRSRC)
lalapps_sned_SOURCES = $(SNEDSRC)
lalapps_corse_SOURCES = $(CORSESRC)
lalapps_skymap_SOURCES = $(SKYMAPSRC)
lalapps_coinj_SOURCES = $(COINJSRC)
lalapps_cbc_stochasticbank_SOURCES = $(STOCHBANKSRC)
endif
endif
endif

if LALFRAME
if LALMETAIO
if LALPULSAR
if LALINSPIRAL
LALPULSAR_PROGS = lalapps_coherent_inspiral
lalapps_coherent_inspiral_SOURCES = $(COHERESRC)
endif
endif
endif
endif

bin_PROGRAMS = \
	$(LALFRAME_LALMETAIO_PROGS) \
	$(LALPULSAR_PROGS)

if HAVE_PYTHON
bin_SCRIPTS = \
	lalapps_inspiral_pipe \
	lalapps_inspiral_hipe \
	lalapps_inspiral_ssipe \
	lalapps_multi_hipe \
	lalapps_trigger_hipe \
	lalapps_inspiral_online_pipe \
	lalapps_inspiral_online_pipe_inject \
	lalapps_ihope_daily_page \
	lalapps_flag_triggers \
	lalapps_check_flag \
	lalapps_check_hipe_times \
	lalapps_ihope \
	lalapps_followup0_pipe \
	lalapps_plot_hipe \
	lalapps_ihope_status \
	lalapps_write_ihope_page \
	lalapps_compare_ini \
	lalapps_link_old_ihope \
	lalapps_generate_upper_limits \
	lalapps_followup0_page \
	lalapps_cbc_pipedown \
	lalapps_cbc_pipedown_ssipe \
	lalapps_glitch_probe \
	lalapps_cohPTF_hipe \
	lalapps_coh_PTF_post_processing \
	lalapps_cbc_sbank \
	lalapps_cbc_sbank_choose_mchirp_boundaries \
	lalapps_cbc_sbank_splitter \
	lalapps_cbc_sbank_merge_sims \
	lalapps_cbc_sbank_pipe \
	lalapps_cbc_sbank_plot_sim \
	lalapps_cbc_sbank_sim

pkgpython_PYTHON = \
	inspiral.py \
	inspiralutils.py
endif

dist_bin_SCRIPTS = \
	lalapps_makeblindinj \
	lalapps_makeblindinj_himass \
	lalapps_cbc_code_compare \
	lalapps_inspiral_ckpttest \
	followup/followup_InspiralDataMover.sh

python_config.sed: Makefile
	@-rm -f python_config.sed
	@echo "s+@PYTHONLIBDIR@+@pkgpythondir@+g" > python_config.sed
	@echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed

ssipe_config.sed: Makefile
	@-rm -f ssipe_config.sed
	@echo "s+@LALAPPSDIR@+@bindir@+g" > ssipe_config.sed

ihope_ssipe.ini: $(srcdir)/ihope_ssipe.ini.in ssipe_config.sed
	$(AM_V_GEN)sed -f ssipe_config.sed $(srcdir)/ihope_ssipe.ini.in > ihope_ssipe.ini

lalapps_inspiral_pipe: $(srcdir)/inspiral_pipe.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/inspiral_pipe.in > lalapps_inspiral_pipe
	@chmod +x lalapps_inspiral_pipe

lalapps_inspiral_hipe: $(srcdir)/inspiral_hipe.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/inspiral_hipe.in > lalapps_inspiral_hipe
	@chmod +x lalapps_inspiral_hipe

lalapps_inspiral_ssipe: $(srcdir)/inspiral_ssipe.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/inspiral_ssipe.in > lalapps_inspiral_ssipe
	@chmod +x lalapps_inspiral_ssipe

lalapps_plot_hipe: $(srcdir)/plot_hipe.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/plot_hipe.in > lalapps_plot_hipe
	@chmod +x lalapps_plot_hipe

lalapps_ihope: $(srcdir)/ihope.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/ihope.in > lalapps_ihope
	@chmod +x lalapps_ihope

lalapps_ihope_status: $(srcdir)/ihope_status.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/ihope_status.in > lalapps_ihope_status
	@chmod +x lalapps_ihope_status

lalapps_write_ihope_page: $(srcdir)/write_ihope_page.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/write_ihope_page.in > lalapps_write_ihope_page
	@chmod +x lalapps_write_ihope_page

lalapps_link_old_ihope: $(srcdir)/link_old_ihope.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/link_old_ihope.in > lalapps_link_old_ihope
	@chmod +x lalapps_link_old_ihope

lalapps_compare_ini: $(srcdir)/compare_ini.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/compare_ini.in > lalapps_compare_ini
	@chmod +x lalapps_compare_ini

lalapps_followup0_pipe: $(srcdir)/followup0_pipe.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/followup0_pipe.in > lalapps_followup0_pipe
	@chmod +x lalapps_followup0_pipe

lalapps_multi_hipe: $(srcdir)/multi_hipe.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/multi_hipe.in > lalapps_multi_hipe
	@chmod +x lalapps_multi_hipe

lalapps_trigger_hipe: $(srcdir)/trigger_hipe.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/trigger_hipe.in > lalapps_trigger_hipe
	@chmod +x lalapps_trigger_hipe

lalapps_cohPTF_hipe: $(srcdir)/cohPTF_hipe.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/cohPTF_hipe.in > lalapps_cohPTF_hipe
	@chmod +x lalapps_cohPTF_hipe

lalapps_check_hipe_times: $(srcdir)/check_hipe_times.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/check_hipe_times.in > lalapps_check_hipe_times
	@chmod +x lalapps_check_hipe_times

lalapps_inspiral_online_pipe: $(srcdir)/inspiral_online_pipe.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/inspiral_online_pipe.in > lalapps_inspiral_online_pipe
	@chmod +x lalapps_inspiral_online_pipe

lalapps_ihope_daily_page: $(srcdir)/ihope_daily_page.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/ihope_daily_page.in > lalapps_ihope_daily_page
	@chmod +x lalapps_ihope_daily_page

lalapps_flag_triggers: $(srcdir)/flag_triggers.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/flag_triggers.in > lalapps_flag_triggers
	@chmod +x lalapps_flag_triggers

lalapps_check_flag: $(srcdir)/check_flag.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/check_flag.in > lalapps_check_flag
	@chmod +x lalapps_check_flag

lalapps_glitch_probe: $(srcdir)/glitch_probe.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/glitch_probe.in > lalapps_glitch_probe
	@chmod +x lalapps_glitch_probe

lalapps_inspiral_online_pipe_inject: $(srcdir)/inspiral_online_pipe_inject.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/inspiral_online_pipe_inject.in > lalapps_inspiral_online_pipe_inject
	@chmod +x lalapps_inspiral_online_pipe_inject

lalapps_generate_upper_limits: $(srcdir)/generate_upper_limits.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/generate_upper_limits.in > lalapps_generate_upper_limits
	@chmod +x lalapps_generate_upper_limits

lalapps_followup0_page: $(srcdir)/followup0_page.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/followup0_page.in > lalapps_followup0_page
	@chmod +x lalapps_followup0_page

lalapps_cbc_pipedown: $(srcdir)/cbc_pipedown.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/cbc_pipedown.in > lalapps_cbc_pipedown
	@chmod +x lalapps_cbc_pipedown

lalapps_cbc_pipedown_ssipe: $(srcdir)/cbc_pipedown_ssipe.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/cbc_pipedown_ssipe.in > lalapps_cbc_pipedown_ssipe
	@chmod +x lalapps_cbc_pipedown_ssipe

lalapps_coh_PTF_post_processing: $(srcdir)/coh_PTF_post_processing.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/coh_PTF_post_processing.in > lalapps_coh_PTF_post_processing
	@chmod +x lalapps_coh_PTF_post_processing

lalapps_cbc_sbank: $(srcdir)/cbc_sbank.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/cbc_sbank.in > lalapps_cbc_sbank
	@chmod +x lalapps_cbc_sbank

lalapps_cbc_sbank_choose_mchirp_boundaries: $(srcdir)/cbc_sbank_choose_mchirp_boundaries.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/cbc_sbank_choose_mchirp_boundaries.in > lalapps_cbc_sbank_choose_mchirp_boundaries
	@chmod +x lalapps_cbc_sbank_choose_mchirp_boundaries

lalapps_cbc_sbank_splitter: $(srcdir)/cbc_sbank_splitter.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/cbc_sbank_splitter.in > lalapps_cbc_sbank_splitter
	@chmod +x lalapps_cbc_sbank_splitter

lalapps_cbc_sbank_merge_sims: $(srcdir)/cbc_sbank_merge_sims.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/cbc_sbank_merge_sims.in > lalapps_cbc_sbank_merge_sims
	@chmod +x lalapps_cbc_sbank_merge_sims

lalapps_cbc_sbank_pipe: $(srcdir)/cbc_sbank_pipe.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/cbc_sbank_pipe.in > lalapps_cbc_sbank_pipe
	@chmod +x lalapps_cbc_sbank_pipe

lalapps_cbc_sbank_plot_sim: $(srcdir)/cbc_sbank_plot_sim.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/cbc_sbank_plot_sim.in > lalapps_cbc_sbank_plot_sim
	@chmod +x lalapps_cbc_sbank_plot_sim

lalapps_cbc_sbank_sim: $(srcdir)/cbc_sbank_sim.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/cbc_sbank_sim.in > lalapps_cbc_sbank_sim
	@chmod +x lalapps_cbc_sbank_sim

pkgdata_DATA = \
	BNSMasses.dat \
	inspsrcs100Mpc.errors \
	ihope_ssipe.ini \
	inspiral_pipe.ini \
	inspiral_hipe.ini \
	online.ini \
	online_inject_bbh.ini \
	cbc_pipedown.ini \
	followup0_pipe.ini \
	plot_hipe.ini \
	write_ihope_page.ini \
	write_ihope_style.css \
	s5_background_H0H1-RDS_R_L1-cbc.txt \
	s5_background_H0H1-RDS_R_L1-seismic-cbc.txt \
	s5_background_H0H2-RDS_R_L1-cbc.txt \
	s5_background_H0H2-RDS_R_L1-seismic-cbc.txt \
	s5_background_H1_hoft_cbc.txt \
	s5_background_H2_hoft_cbc.txt \
	s5_background_L0L1-RDS_R_L1-cbc.txt \
	s5_background_L0L1-RDS_R_L1-seismic-cbc.txt \
	s5_background_L1_hoft_cbc.txt \
	s5_background_V1-raw-cbc.txt \
	s5_background_V1-raw-seismic-cbc.txt \
	s5_background_V1_hoft_cbc.txt \
	s5_foreground_H0H1-RDS_R_L1-cbc.txt \
	s5_foreground_H0H1-RDS_R_L1-seismic-cbc.txt \
	s5_foreground_H0H2-RDS_R_L1-cbc.txt \
	s5_foreground_H0H2-RDS_R_L1-seismic-cbc.txt \
	s5_foreground_H1_hoft_cbc.txt \
	s5_foreground_H2_hoft_cbc.txt \
	s5_foreground_L0L1-RDS_R_L1-cbc.txt \
	s5_foreground_L0L1-RDS_R_L1-seismic-cbc.txt \
	s5_foreground_L1_hoft_cbc.txt \
	s5_foreground_V1-raw-cbc.txt \
	s5_foreground_V1-raw-seismic-cbc.txt \
	s5_foreground_V1_hoft_cbc.txt \
	s6_background_H0H1-RDS_R_L1-cbc.txt \
	s6_background_H0H1-RDS_R_L1-seismic-cbc.txt \
	s6_background_H1_hoft_cbc.txt \
	s6_background_L0L1-RDS_R_L1-cbc.txt \
	s6_background_L0L1-RDS_R_L1-seismic-cbc.txt \
	s6_background_L1_hoft_cbc.txt \
	s6_background_V1-raw-cbc.txt \
	s6_background_V1-raw-seismic-cbc.txt \
	s6_background_V1_hoft_cbc.txt \
	s6_foreground_H0H1-RDS_R_L1-cbc.txt \
	s6_foreground_H0H1-RDS_R_L1-seismic-cbc.txt \
	s6_foreground_H1_hoft_cbc.txt \
	s6_foreground_L0L1-RDS_R_L1-cbc.txt \
	s6_foreground_L0L1-RDS_R_L1-seismic-cbc.txt \
	s6_foreground_L1_hoft_cbc.txt \
	s6_foreground_V1-raw-cbc.txt \
	s6_foreground_V1-raw-seismic-cbc.txt \
	s6_foreground_V1_hoft_cbc.txt \
	ratioTest.pickle \
	remove_h1h2.sql \
	cluster.sql \
	simplify.sql \
	mcmcsummary.R

EXTRA_DIST = \
	$(INSPIRALSRC) \
	$(TMPLTBANKSRC) \
	$(SIRESRC) \
	$(COIRESRC) \
	$(SIVASRC) \
	$(MINJSRC) \
	$(INCASRC) \
	$(COHERESRC) \
	$(INSPINJSRC) \
	$(SPLITBANKSRC) \
	inspiral.py \
	inspiralutils.py \
	inspiral_pipe.in \
	inspiral_online_pipe.in \
	inspiral_hipe.in \
	inspiral_ssipe.in \
	cohPTF_hipe.in \
	multi_hipe.in \
	trigger_hipe.in \
	ihope.in \
	ihope_ssipe.ini.in \
	ihope_status.in \
	inspiral_online_pipe_inject.in \
	check_hipe_times.in \
	followup0_pipe.in \
	plot_hipe.in \
	write_ihope_page.in \
	compare_ini.in \
	link_old_ihope.in \
	generate_upper_limits.in \
	followup0_page.in \
	cbc_pipedown.in \
	cbc_pipedown_ssipe.in \
	ihope_daily_page.in \
	check_flag.in \
	flag_triggers.in \
	glitch_probe.in \
	coh_PTF_post_processing.in \
	cbc_sbank.in \
	cbc_sbank_choose_mchirp_boundaries.in \
	cbc_sbank_merge_sims.in \
	cbc_sbank_pipe.in \
	cbc_sbank_plot_sim.in \
	cbc_sbank_sim.in \
	cbc_sbank_splitter.in \
	$(DATA)

AM_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_builddir)/src/lalapps -DPKG_DATA_DIR='"$(pkgdatadir)/"'
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

CLEANFILES = \
	$(bin_SCRIPTS) \
	python_config.sed \
	ssipe_config.sed \
	ihope_ssipe.ini

MOSTLYCLEANFILES = \
	*INSPIRAL*.xml \
	*INSPIRAL*.gwf \
	*TMPLTBANK*.gwf \
	*TMPLTBANK*.xml

BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
