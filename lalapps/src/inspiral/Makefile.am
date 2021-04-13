BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

TMPLTBANKSRC = tmpltbank.c inspiralutils.c inspiral.h
INSPINJSRC = inspinj.c inspiralutils.c inspiral.h
BLINDINJSRC = blindinj.c
SPLITBANKSRC = splitbank.c
TRIGSCANSRC = trigscan.c
INSPFRINJSRC = inspfrinj.c
GWF2XMLSRC = gwf2xml.c
SPININJSRC = spininj.c
RANDOMBANKSRC = randombank.c
NINJASRC = ninja.c
FRNINJASRC = fr_ninja.c
CALCEXPSNRSRC = calcexpsnr.c
COINJSRC = coinj.c
STOCHBANKSRC = stochasticbank.c

lalapps_tmpltbank_SOURCES = $(TMPLTBANKSRC)
lalapps_inspinj_SOURCES = $(INSPINJSRC)
lalapps_blindinj_SOURCES = $(BLINDINJSRC)
lalapps_splitbank_SOURCES = $(SPLITBANKSRC)
lalapps_trigscan_SOURCES = $(TRIGSCANSRC)
lalapps_spininj_SOURCES = $(SPININJSRC)
lalapps_inspfrinj_SOURCES = $(INSPFRINJSRC)
lalapps_gwf2xml_SOURCES = $(GWF2XMLSRC)
lalapps_randombank_SOURCES = $(RANDOMBANKSRC)
lalapps_ninja_SOURCES = $(NINJASRC)
lalapps_fr_ninja_SOURCES = $(FRNINJASRC)
lalapps_calcexpsnr_SOURCES = $(CALCEXPSNRSRC)
lalapps_coinj_SOURCES = $(COINJSRC)
lalapps_cbc_stochasticbank_SOURCES = $(STOCHBANKSRC)

bin_PROGRAMS = \
	lalapps_splitbank \
	lalapps_blindinj \
	lalapps_trigscan \
	lalapps_spininj \
	lalapps_randombank \
	lalapps_fr_ninja \
	lalapps_coinj \
	lalapps_cbc_stochasticbank \
	$(END_OF_LIST)

if FRAMEL
bin_PROGRAMS += \
	lalapps_tmpltbank \
	lalapps_inspinj \
	lalapps_inspfrinj \
	lalapps_gwf2xml \
	lalapps_ninja \
	lalapps_calcexpsnr \
	$(END_OF_LIST)
endif

if HAVE_PYTHON
pybin_scripts = \
	lalapps_inspinjfind \
	lalapps_thinca \
	lalapps_cbc_sbank \
	lalapps_cbc_sbank_choose_mchirp_boundaries \
	lalapps_cbc_sbank_hdf5_choose_mchirp_boundaries \
	lalapps_cbc_sbank_hdf5_bankcombiner \
	lalapps_cbc_sbank_pipe \
	lalapps_cbc_sbank_plot_sim \
	lalapps_cbc_sbank_sim \
	lalapps_cbc_sbank_sim_pipe \
	lalapps_make_nr_hdf_catalog \
	lalapps_run_sqlite \
	$(END_OF_LIST)

pkgpython_PYTHON = \
	inspiral.py \
	inspiralutils.py \
	$(END_OF_LIST)
endif

dist_bin_SCRIPTS = \
	lalapps_makeblindinj \
	lalapps_makeblindinj_himass \
	$(END_OF_LIST)

pkgdata_DATA = \
	BNSMasses.dat \
	inspsrcs100Mpc.errors \
	$(END_OF_LIST)

MOSTLYCLEANFILES += \
	*INSPIRAL*.xml \
	*INSPIRAL*.gwf \
	*TMPLTBANK*.gwf \
	*TMPLTBANK*.xml \
	$(END_OF_LIST)

EXTRA_DIST += \
	$(DATA) \
	$(INSPINJSRC) \
	$(SPLITBANKSRC) \
	$(TMPLTBANKSRC) \
	effective_distance_comparison.png \
	inspiral_programs.dox \
	$(END_OF_LIST)
