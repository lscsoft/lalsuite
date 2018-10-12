include $(top_srcdir)/gnuscripts/lalapps.am

MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

RINGPROG = \
	lalapps_coincringread \
	lalapps_coh_PTF_inspiral \
	lalapps_coh_PTF_spin_checker

bin_PROGRAMS = $(RINGPROG)

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
lalapps_coh_PTF_inspiral_SOURCES = \
	coh_PTF_inspiral.c \
	$(COH_PTF_SOURCES) \
	$(GENERAL_SOURCES)
lalapps_coh_PTF_spin_checker_SOURCES = \
	coh_PTF_spin_checker.c \
	$(COH_PTF_SOURCES) \
	$(GENERAL_SOURCES)
lalapps_coincringread_SOURCES = coincringread.c

if HAVE_PYTHON
pkgpython_PYTHON = ring.py
pybin_scripts = \
	lalapps_ring_pipe \
	lalapps_ringcorse
endif

pkgdata_DATA = ring_test.ini ringpost_remove_cluster.sql
EXTRA_DIST += $(DATA)
