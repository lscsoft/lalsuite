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

noinst_PROGRAMS = \
	genSFTWindows \
	$(END_OF_LIST)

EXTRA_PROGRAMS =

lalpulsar_Makefakedata_v4_SOURCES = \
	makefakedata_v4.c \
	$(END_OF_LIST)

lalpulsar_Makefakedata_v5_SOURCES = \
	makefakedata_v5.c \
	$(END_OF_LIST)

genSFTWindows_SOURCES = \
	SFTWindowFuncs.c \
	genSFTWindows.c \
	$(END_OF_LIST)

if LALFRAME

bin_PROGRAMS += \
	lalpulsar_MakeSFTs \
	$(END_OF_LIST)

lalpulsar_MakeSFTs_SOURCES = \
	MakeSFTs.c \
	$(END_OF_LIST)

EXTRA_PROGRAMS += \
	lalpulsar_sw_inj_frames \
	$(END_OF_LIST)

lalpulsar_sw_inj_frames_SOURCES = \
	sw_inj_frames.c \
	$(END_OF_LIST)

endif

if HAVE_PYTHON

pybin_scripts = \
	lalpulsar_CopyPublicSFTs \
	lalpulsar_MakeSFTDAG \
	$(END_OF_LIST)

endif

# Add shell test scripts to this variable
test_scripts += testMakeSFTs.sh
test_scripts += testMakefakedata_v4.sh
test_scripts += testMakefakedata_v4_hwinjection.sh
test_scripts += testMakefakedata_v5.sh
if HAVE_PYTHON
test_scripts += testCopyPublicSFTs.sh
test_scripts += testMakeSFTDAG.sh
test_scripts += testSFTWindows.py
endif

# Add any helper programs required by tests to this variable
test_helpers += compareTS
test_helpers += genSFTWindows

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
