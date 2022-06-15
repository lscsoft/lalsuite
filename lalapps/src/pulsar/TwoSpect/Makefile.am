BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps_pulsar.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

noinst_PROGRAMS = \
	lalapps_TwoSpect \
	lalapps_TwoSpectTemplateBank \
	$(END_OF_LIST)

EXTRA_PROGRAMS = \
	compareCandidates \
	computeSignalDetector \
	skygridsetup \
	testVectorMath \
	$(END_OF_LIST)

lalapps_TwoSpect_SOURCES = \
	IHS.c\
	IHS.h \
	SFTfunctions.c \
	SFTfunctions.h \
	TwoSpect.c \
	TwoSpect.h \
	TwoSpectSpecFunc.c \
	TwoSpectSpecFunc.h \
	TwoSpectTypes.h \
	antenna.c \
	antenna.h \
	candidates.c \
	candidates.h \
	cdfdist.c \
	cdfdist.h \
	cdfwchisq.c \
	cdfwchisq.h \
	falsealarm.c \
	falsealarm.h \
	statistics.c \
	statistics.h \
	templates.c \
	templates.h \
	upperlimits.c \
	upperlimits.h \
	vectormath.c \
	vectormath.h \
	$(END_OF_LIST)

lalapps_TwoSpectTemplateBank_SOURCES = \
	TwoSpectSpecFunc.c \
	TwoSpectSpecFunc.h \
	TwoSpectTemplateBank.c \
	statistics.c \
	statistics.h \
	templates.c \
	templates.h \
	vectormath.c \
	vectormath.h \
	$(END_OF_LIST)

skygridsetup_SOURCES = \
	antenna.c \
	antenna.h \
	helperprograms/skygridsetup.c \
	$(END_OF_LIST)

compareCandidates_SOURCES = \
	helperprograms/compareCandidates.c \
	$(END_OF_LIST)

testVectorMath_SOURCES = \
	helperprograms/testVectorMath.c \
	vectormath.c \
	vectormath.h \
	$(END_OF_LIST)

computeSignalDetector_SOURCES = \
	TwoSpectSpecFunc.c \
	TwoSpectSpecFunc.h \
	helperprograms/computeSignalDetector.c \
	$(END_OF_LIST)

# Add shell test scripts to this variable
test_scripts +=

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
