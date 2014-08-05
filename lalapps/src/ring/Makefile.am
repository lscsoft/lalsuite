BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalapps.am

if FRAME
if LALFRAME
if LALMETAIO
if LALINSPIRAL
RINGPROG = \
	lalapps_ring \
	lalapps_rinj \
	lalapps_ringread \
	lalapps_rinca \
	lalapps_coincringread \
	lalapps_collatesim \
	lalapps_coh_PTF_inspiral \
	lalapps_coh_PTF_spin_checker
endif
endif
endif
endif

bin_PROGRAMS = $(RINGPROG)

if FRAME
if LALFRAME
if LALMETAIO
if LALINSPIRAL
GENERAL_SOURCES = \
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
	segment.h
COH_PTF_SOURCES = \
	coh_PTF.h \
	coh_PTF_option.c \
	coh_PTF_template.c \
	coh_PTF_output.c \
	coh_PTF_bankveto.c \
	coh_PTF_utils.c
lalapps_ring_SOURCES = \
	ring.c \
	ring.h \
	ring_option.c \
	ring_filter.c \
	ring_output.c \
	$(GENERAL_SOURCES)
lalapps_coh_PTF_inspiral_SOURCES = \
	coh_PTF_inspiral.c \
	$(COH_PTF_SOURCES) \
	$(GENERAL_SOURCES)
lalapps_coh_PTF_spin_checker_SOURCES = \
	coh_PTF_spin_checker.c \
	$(COH_PTF_SOURCES) \
	$(GENERAL_SOURCES)
lalapps_rinj_SOURCES = rinj.c
lalapps_ringread_SOURCES = ringread.c
lalapps_rinca_SOURCES = rinca.c
lalapps_coincringread_SOURCES = coincringread.c
lalapps_collatesim_SOURCES = collatesim.c
endif
endif
endif
endif

if HAVE_PYTHON
pkgpython_PYTHON = ring.py
pybin_SCRIPTS = \
	lalapps_ring_pipe \
	lalapps_ringcorse
endif

pkgdata_DATA = ring_test.ini ringpost_remove_cluster.sql
EXTRA_DIST += $(DATA)
