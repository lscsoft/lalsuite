BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps_pulsar.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalapps_DriveHoughMulti \
	$(END_OF_LIST)

noinst_PROGRAMS = \
	lalapps_HoughValidate \
	lalapps_HoughValidateAM \
	lalapps_MCInjectHoughMulti \
	lalapps_MCInjectHoughMultiChi2Test \
	lalapps_MultiWeights \
	lalapps_ValidateChi2Test \
	lalapps_ValidateHoughMulti \
	$(END_OF_LIST)

lalapps_HoughValidate_SOURCES = \
	DriveHoughColor.h \
	HoughValidate.c \
	MCInjectComputeHough.h \
	PeakSelect.c \
	PeakSelect.h \
	$(END_OF_LIST)

lalapps_DriveHoughMulti_SOURCES = \
	DriveHoughColor.h \
	DriveHoughMulti.c \
	FstatToplist.c \
	FstatToplist.h \
	PeakSelect.c \
	PeakSelect.h \
	$(END_OF_LIST)

lalapps_ValidateHoughMulti_SOURCES = \
	DriveHoughColor.h \
	MCInjectHoughMulti.h \
	PeakSelect.c \
	PeakSelect.h \
	ValidateHoughMulti.c \
	$(END_OF_LIST)

lalapps_MultiWeights_SOURCES = \
	DriveHoughColor.h \
	MultiWeights.c \
	$(END_OF_LIST)

lalapps_MCInjectHoughMulti_SOURCES = \
	DriveHoughColor.h \
	MCInjectHoughMulti.c \
	MCInjectHoughMulti.h \
	PeakSelect.c \
	PeakSelect.h \
	$(END_OF_LIST)

lalapps_ValidateChi2Test_SOURCES = \
	DriveHoughColor.h \
	MCInjectHoughMulti.h \
	PeakSelect.c \
	PeakSelect.h \
	ValidateHoughMultiChi2Test.c \
	$(END_OF_LIST)

lalapps_MCInjectHoughMultiChi2Test_SOURCES = \
	DriveHoughColor.h \
	MCInjectHoughMulti.h \
	MCInjectHoughMultiChi2Test.c \
	PeakSelect.c \
	PeakSelect.h \
	$(END_OF_LIST)

lalapps_HoughValidateAM_SOURCES = \
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
