BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

TMPLTBANKSRC = tmpltbank.c inspiralutils.c inspiral.h
INSPINJSRC = inspinj.c inspiralutils.c inspiral.h LALExtTriggerTableFromLIGOLw.h LALExtTriggerTableFromLIGOLw.c
BLINDINJSRC = blindinj.c
SPLITBANKSRC = splitbank.c
INSPFRINJSRC = inspfrinj.c
GWF2XMLSRC = gwf2xml.c
SPININJSRC = spininj.c
RANDOMBANKSRC = randombank.c
NINJASRC = ninja.c
FRNINJASRC = fr_ninja.c
COINJSRC = coinj.c
STOCHBANKSRC = stochasticbank.c

lalapps_tmpltbank_SOURCES = $(TMPLTBANKSRC)
lalapps_inspinj_SOURCES = $(INSPINJSRC)
lalapps_blindinj_SOURCES = $(BLINDINJSRC)
lalapps_splitbank_SOURCES = $(SPLITBANKSRC)
lalapps_spininj_SOURCES = $(SPININJSRC)
lalapps_inspfrinj_SOURCES = $(INSPFRINJSRC)
lalapps_gwf2xml_SOURCES = $(GWF2XMLSRC)
lalapps_randombank_SOURCES = $(RANDOMBANKSRC)
lalapps_ninja_SOURCES = $(NINJASRC)
lalapps_fr_ninja_SOURCES = $(FRNINJASRC)
lalapps_coinj_SOURCES = $(COINJSRC)
lalapps_cbc_stochasticbank_SOURCES = $(STOCHBANKSRC)

bin_PROGRAMS = \
	lalapps_splitbank \
	lalapps_blindinj \
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
	$(END_OF_LIST)
endif

if HAVE_PYTHON
pybin_scripts = \
	$(END_OF_LIST)

pkgpython_PYTHON = \
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
