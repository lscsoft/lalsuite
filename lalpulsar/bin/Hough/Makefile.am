BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
man1_MANS =
include $(top_srcdir)/gnuscripts/lalpulsar_bin.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalpulsar_DriveHoughMulti \
	$(END_OF_LIST)

noinst_PROGRAMS = \
	lalpulsar_HoughValidate \
	lalpulsar_HoughValidateAM \
	lalpulsar_MCInjectHoughMulti \
	lalpulsar_MCInjectHoughMultiChi2Test \
	lalpulsar_MultiWeights \
	lalpulsar_ValidateChi2Test \
	lalpulsar_ValidateHoughMulti \
	$(END_OF_LIST)

lalpulsar_HoughValidate_SOURCES = \
	DriveHoughColor.h \
	HoughValidate.c \
	MCInjectComputeHough.h \
	PeakSelect.c \
	PeakSelect.h \
	$(END_OF_LIST)

lalpulsar_DriveHoughMulti_SOURCES = \
	DriveHoughColor.h \
	DriveHoughMulti.c \
	FstatToplist.c \
	FstatToplist.h \
	PeakSelect.c \
	PeakSelect.h \
	$(END_OF_LIST)

lalpulsar_ValidateHoughMulti_SOURCES = \
	DriveHoughColor.h \
	MCInjectHoughMulti.h \
	PeakSelect.c \
	PeakSelect.h \
	ValidateHoughMulti.c \
	$(END_OF_LIST)

lalpulsar_MultiWeights_SOURCES = \
	DriveHoughColor.h \
	MultiWeights.c \
	$(END_OF_LIST)

lalpulsar_MCInjectHoughMulti_SOURCES = \
	DriveHoughColor.h \
	MCInjectHoughMulti.c \
	MCInjectHoughMulti.h \
	PeakSelect.c \
	PeakSelect.h \
	$(END_OF_LIST)

lalpulsar_ValidateChi2Test_SOURCES = \
	DriveHoughColor.h \
	MCInjectHoughMulti.h \
	PeakSelect.c \
	PeakSelect.h \
	ValidateHoughMultiChi2Test.c \
	$(END_OF_LIST)

lalpulsar_MCInjectHoughMultiChi2Test_SOURCES = \
	DriveHoughColor.h \
	MCInjectHoughMulti.h \
	MCInjectHoughMultiChi2Test.c \
	PeakSelect.c \
	PeakSelect.h \
	$(END_OF_LIST)

lalpulsar_HoughValidateAM_SOURCES = \
	DriveHoughColor.h \
	HoughValidateAM.c \
	MCInjectHoughS2.h \
	PeakSelect.c \
	PeakSelect.h \
	SFTbin.c \
	SFTbin.h \
	$(END_OF_LIST)

# Add shell test scripts to this variable
if HAVE_PYTHON
test_scripts += testDriveHoughMulti.py
endif

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
