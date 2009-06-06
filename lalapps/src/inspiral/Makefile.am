## Process this file with automake to produce Makefile.in
# $Id$

SUBDIRS = posterior

INSPIRALSRC 	= inspiral.c inspiralutils.c inspiral.h
TMPLTBANKSRC 	= tmpltbank.c inspiralutils.c inspiral.h
COHERESRC 	= coherent_inspiral.c
COHIRESRC	= cohire.c
COHEREBANK 	= coherentbank.c
INCASRC 	= inca.c
THINCASRC 	= thinca.c
SIRESRC 	= sire.c
COIRESRC 	= coire.c
SIVASRC 	= siva.c
INSPINJSRC 	= inspinj.c
INJCUTSRC 	= injcut.c
MINJSRC 	= minj.c
BBHINJSRC 	= bbhinj.c
BLINDINJSRC 	= blindinj.c
SPLITBANKSRC 	= splitbank.c
TRIGBANKSRC 	= trigbank.c
INSPFRINJSRC 	= inspfrinj.c
GWF2XMLSRC 	= gwf2xml.c
SPININJSRC 	= spininj.c
COINEXTSRC 	= coinext.c
RANDOMBANKSRC 	= randombank.c
NINJASRC	= ninja.c
FRNINJASRC	= fr_ninja.c
MDCNINJASRC	= mdc_ninja.c inspiralutils.c inspiral.h
CALCEXPSNRSRC	= calcexpsnr.c
SNEDSRC		= sned.c
CORSESRC	= corse.c
SKYMAPSRC	= skymap.c
COINJSRC	= coinj.c

## make sure LALAPPSgitID.h is re-generated at every 'make'
.PHONY: gitID
gitID:
	cd $(top_builddir)/misc && $(MAKE) gitID
BUILT_SOURCES = gitID

if LALFRAME
if LALMETAIO
bin_PROGRAMS = lalapps_inspiral lalapps_tmpltbank \
               lalapps_inca lalapps_sire lalapps_coire lalapps_siva \
               lalapps_inspinj lalapps_injcut lalapps_splitbank \
               lalapps_minj lalapps_bbhinj lalapps_blindinj \
               lalapps_inspfrinj lalapps_gwf2xml \
               lalapps_thinca lalapps_trigbank \
               lalapps_spininj lalapps_coherent_inspiral lalapps_cohire \
               lalapps_coherentbank lalapps_coinext \
               lalapps_randombank lalapps_ninja lalapps_fr_ninja \
               lalapps_mdc_ninja lalapps_calcexpsnr lalapps_sned \
               lalapps_corse lalapps_skymap lalapps_coinj

lalapps_inspiral_SOURCES 	= $(INSPIRALSRC)
lalapps_coherent_inspiral_SOURCES = $(COHERESRC)
lalapps_cohire_SOURCES          = $(COHIRESRC)
lalapps_coherentbank_SOURCES 	= $(COHEREBANK)
lalapps_tmpltbank_SOURCES 	= $(TMPLTBANKSRC)
lalapps_inca_SOURCES 		= $(INCASRC)
lalapps_thinca_SOURCES 		= $(THINCASRC)
lalapps_sire_SOURCES 		= $(SIRESRC)
lalapps_coire_SOURCES 		= $(COIRESRC)
lalapps_siva_SOURCES 		= $(SIVASRC)
lalapps_inspinj_SOURCES 	= $(INSPINJSRC)
lalapps_injcut_SOURCES 		= $(INJCUTSRC)
lalapps_minj_SOURCES 		= $(MINJSRC)
lalapps_bbhinj_SOURCES 		= $(BBHINJSRC)
lalapps_blindinj_SOURCES 	= $(BLINDINJSRC)
lalapps_splitbank_SOURCES 	= $(SPLITBANKSRC)
lalapps_trigbank_SOURCES 	= $(TRIGBANKSRC)
lalapps_spininj_SOURCES 	= $(SPININJSRC)
lalapps_inspfrinj_SOURCES 	= $(INSPFRINJSRC)
lalapps_gwf2xml_SOURCES 	= $(GWF2XMLSRC)
lalapps_coinext_SOURCES         = $(COINEXTSRC)
lalapps_randombank_SOURCES 	= $(RANDOMBANKSRC)
lalapps_ninja_SOURCES		= $(NINJASRC)
lalapps_fr_ninja_SOURCES	= $(FRNINJASRC)
lalapps_mdc_ninja_SOURCES	= $(MDCNINJASRC)
lalapps_calcexpsnr_SOURCES	= $(CALCEXPSNRSRC)
lalapps_sned_SOURCES		= $(SNEDSRC)
lalapps_corse_SOURCES		= $(CORSESRC)
lalapps_skymap_SOURCES		= $(SKYMAPSRC)
lalapps_coinj_SOURCES		= $(COINJSRC)
else
bin_PROGRAMS = 
endif
else
bin_PROGRAMS = 
endif

if HAVE_PYTHON
bin_SCRIPTS = lalapps_inspiral_pipe lalapps_inspiral_hipe \
	      lalapps_multi_hipe lalapps_trigger_hipe\
              lalapps_inspiral_online_pipe \
              lalapps_inspiral_online_pipe_inject \
              lalapps_check_hipe_times lalapps_ihope \
              lalapps_followup_pipe lalapps_plot_hipe \
              lalapps_ihope_status lalapps_write_ihope_page \
	      lalapps_compare_ini lalapps_link_old_ihope \
	      lalapps_generate_upper_limits lalapps_followup_page \
	      lalapps_makeblindinj lalapps_makeblindinj_himass

pkgpython_DATA = inspiral.py inspiralutils.py
else
bin_SCRIPTS = lalapps_makeblindinj lalapps_makeblindinj_himass
endif

lalapps_makeblindinj:
	chmod a+x lalapps_makeblindinj

lalapps_makeblindinj_himass:
	chmod a+x lalapps_makeblindinj_himass

python_config.sed:
	-rm -f python_config.sed
	echo "s+@PYTHONLIBDIR@+@pkgpythondir@+g" > python_config.sed
	echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed

lalapps_inspiral_pipe: $(srcdir)/inspiral_pipe.in python_config.sed
	sed -f python_config.sed $(srcdir)/inspiral_pipe.in > lalapps_inspiral_pipe 
	chmod +x lalapps_inspiral_pipe
lalapps_inspiral_hipe: $(srcdir)/inspiral_hipe.in python_config.sed
	sed -f python_config.sed $(srcdir)/inspiral_hipe.in > lalapps_inspiral_hipe
	chmod +x lalapps_inspiral_hipe
lalapps_plot_hipe: $(srcdir)/plot_hipe.in python_config.sed
	sed -f python_config.sed $(srcdir)/plot_hipe.in > lalapps_plot_hipe
	chmod +x lalapps_plot_hipe
lalapps_ihope: $(srcdir)/ihope.in python_config.sed
	sed -f python_config.sed $(srcdir)/ihope.in > lalapps_ihope
	chmod +x lalapps_ihope
lalapps_ihope_status: $(srcdir)/ihope_status.in python_config.sed
	sed -f python_config.sed $(srcdir)/ihope_status.in > lalapps_ihope_status
	chmod +x lalapps_ihope_status
lalapps_write_ihope_page: $(srcdir)/write_ihope_page.in python_config.sed
	sed -f python_config.sed $(srcdir)/write_ihope_page.in > lalapps_write_ihope_page
	chmod +x lalapps_write_ihope_page
lalapps_link_old_ihope: $(srcdir)/link_old_ihope.in python_config.sed
	sed -f python_config.sed $(srcdir)/link_old_ihope.in > lalapps_link_old_ihope
	chmod +x lalapps_link_old_ihope
lalapps_compare_ini: $(srcdir)/compare_ini.in python_config.sed
	sed -f python_config.sed $(srcdir)/compare_ini.in > lalapps_compare_ini
	chmod +x lalapps_compare_ini
lalapps_followup_pipe: $(srcdir)/followup_pipe.in python_config.sed
lalapps_followup_pipe: $(srcdir)/followup_pipe.in python_config.sed
	sed -f python_config.sed $(srcdir)/followup_pipe.in > lalapps_followup_pipe
	chmod +x lalapps_followup_pipe
lalapps_multi_hipe: $(srcdir)/multi_hipe.in python_config.sed
	sed -f python_config.sed $(srcdir)/multi_hipe.in > lalapps_multi_hipe
	chmod +x lalapps_multi_hipe
lalapps_trigger_hipe: $(srcdir)/trigger_hipe.in python_config.sed
	sed -f python_config.sed $(srcdir)/trigger_hipe.in > lalapps_trigger_hipe
	chmod +x lalapps_trigger_hipe
lalapps_check_hipe_times: $(srcdir)/check_hipe_times.in python_config.sed
	sed -f python_config.sed $(srcdir)/check_hipe_times.in > lalapps_check_hipe_times
	chmod +x lalapps_check_hipe_times
lalapps_inspiral_online_pipe: $(srcdir)/inspiral_online_pipe.in python_config.sed
	sed -f python_config.sed $(srcdir)/inspiral_online_pipe.in > lalapps_inspiral_online_pipe
	chmod +x lalapps_inspiral_online_pipe
lalapps_inspiral_online_pipe_inject: $(srcdir)/inspiral_online_pipe_inject.in python_config.sed
	sed -f python_config.sed $(srcdir)/inspiral_online_pipe_inject.in > lalapps_inspiral_online_pipe_inject
	chmod +x lalapps_inspiral_online_pipe_inject
lalapps_generate_upper_limits: $(srcdir)/generate_upper_limits.in python_config.sed
	sed -f python_config.sed $(srcdir)/generate_upper_limits.in > lalapps_generate_upper_limits
	chmod +x lalapps_generate_upper_limits
lalapps_followup_page: $(srcdir)/followup_page.in python_config.sed
	sed -f python_config.sed $(srcdir)/followup_page.in > lalapps_followup_page
	chmod +x lalapps_followup_page


pkgdata_DATA = BNSMasses.dat inspsrcs100Mpc.complete inspsrcs100Mpc.errors inspiral_pipe.ini inspiral_hipe.ini online.ini online_inject_bbh.ini ihope.ini followup_pipe.ini plot_hipe.ini write_ihope_page.ini write_ihope_style.css
EXTRA_DIST = $(INSPIRALSRC) $(TMPLTBANKSRC) \
             $(SIRESRC) $(COIRESRC) $(SIVASRC) $(MINJSRC) \
             $(INCASRC) $(COHERESRC) $(INSPINJSRC) \
             $(SPLITBANKSRC) inspiral.py inspiralutils.py inspiral_pipe.in \
             inspiral_online_pipe.in inspiral_hipe.in \
	     multi_hipe.in trigger_hipe.in ihope.in ihope_status.in\
             inspiral_online_pipe_inject.in check_hipe_times.in \
	     followup_pipe.in plot_hipe.in  write_ihope_page.in\
	     compare_ini.in link_old_ihope.in generate_upper_limits.in\
	     followup_page.in lalapps_makeblindinj lalapps_makeblindinj_himass\
             $(DATA)

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_srcdir)/misc
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

CLEANFILES = $(bin_SCRIPTS) python_config.sed

MOSTLYCLEANFILES = *INSPIRAL*.xml *INSPIRAL*.gwf \
                   *TMPLTBANK*.gwf *TMPLTBANK*.xml
