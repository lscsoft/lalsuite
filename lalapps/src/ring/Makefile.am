BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

bin_PROGRAMS = \
	lalapps_coh_PTF_inspiral \
	lalapps_coh_PTF_spin_checker

COH_PTF_SOURCES = \
	errutil.c \
	errutil.h \
	gpstime.c \
	gpstime.h \
	getdata.c \
	getdata.h \
	injsgnl.c \
	injsgnl.h \
	getresp.c \
	getresp.h \
	spectrm.c \
	spectrm.h \
	segment.c \
	segment.h \
	coh_PTF.h \
	coh_PTF_option.c \
	coh_PTF_template.c \
	coh_PTF_output.c \
	coh_PTF_bankveto.c \
	coh_PTF_utils.c \
	FindChirpDatatypes.h \
	FindChirpIMRSimulation.h \
	FindChirpIMRSimulation.c \
	FindChirpTD.h \
	FindChirpTDTemplate.c \
	FindChirpSP.h \
	FindChirpSPTemplate.c \
	GenerateRing.h \
	GenerateRing.c \
	Inject.h \
	InjectTimeSeries.c
lalapps_coh_PTF_inspiral_SOURCES = \
	coh_PTF_inspiral.c \
	$(COH_PTF_SOURCES)
lalapps_coh_PTF_spin_checker_SOURCES = \
	coh_PTF_spin_checker.c \
	$(COH_PTF_SOURCES)
