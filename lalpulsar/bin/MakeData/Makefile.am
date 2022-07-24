BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalpulsar_bin.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalpulsar_Makefakedata_v4 \
	lalpulsar_Makefakedata_v5 \
	$(END_OF_LIST)

EXTRA_PROGRAMS =

lalpulsar_Makefakedata_v4_SOURCES = \
	makefakedata_v4.c \
	$(END_OF_LIST)

lalpulsar_Makefakedata_v5_SOURCES = \
	makefakedata_v5.c \
	$(END_OF_LIST)

if LALFRAME

bin_PROGRAMS += \
	lalpulsar_MakeSFTs \
	$(END_OF_LIST)

lalpulsar_MakeSFTs_SOURCES = \
	MakeSFTs.c \
	SFTWindowFuncs.c \
	$(END_OF_LIST)

lalpulsar_MakeSFTs_CPPFLAGS = $(AM_CPPFLAGS)

if PSS

lalpulsar_MakeSFTs_SOURCES += XLALPSSInterface.c
lalpulsar_MakeSFTs_CPPFLAGS += -DPSS_ENABLED

endif

EXTRA_PROGRAMS += \
	lalpulsar_sw_inj_frames \
	$(END_OF_LIST)

lalpulsar_sw_inj_frames_SOURCES = \
	sw_inj_frames.c \
	$(END_OF_LIST)

endif

dist_bin_SCRIPTS = \
	lalpulsar_CopySFTs \
	$(END_OF_LIST)

if HAVE_PYTHON

pybin_scripts = \
	lalpulsar_MakeSFTDAG \
	$(END_OF_LIST)

endif

# Add shell test scripts to this variable
test_scripts += testMakeSFTs.sh
test_scripts += testMakefakedata_v4.sh
test_scripts += testMakefakedata_v4_hwinjection.sh
test_scripts += testMakefakedata_v5.sh
if HAVE_PYTHON
test_scripts += testMakeSFTDAG.sh
endif

# Add any helper programs required by tests to this variable
test_helpers += compareTS

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
