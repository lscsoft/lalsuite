BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/lalapps.am

bin_PROGRAMS = lalapps_TwoSpect

lalapps_TwoSpect_SOURCES = statistics.c fastchisqinv.c IHS.c candidates.c templates.c antenna.c \
                           vectormath.c cdfwchisq.c upperlimits.c TwoSpect.c \
			   TwoSpectTypes.h statistics.h fastchisqinv.h IHS.h candidates.h templates.h \
                           antenna.h vectormath.h cdfwchisq.h upperlimits.h TwoSpect.h

EXTRA_PROGRAMS = skygridsetup compareCandidates testVectorMath

skygridsetup_SOURCES = helperprograms/skygridsetup.c antenna.c antenna.h

compareCandidates_SOURCES = helperprograms/compareCandidates.c

testVectorMath_SOURCES = helperprograms/testVectorMath.c vectormath.c vectormath.h
