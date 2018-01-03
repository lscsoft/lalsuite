include $(top_srcdir)/gnuscripts/lalapps.am

bin_PROGRAMS = lalapps_TwoSpect lalapps_TwoSpectTemplateBank

lalapps_TwoSpect_SOURCES = statistics.c falsealarm.c TwoSpectSpecFunc.c IHS.c candidates.c \
                           templates.c antenna.c vectormath.c cdfwchisq.c upperlimits.c cdfdist.c \
                           SFTfunctions.c TwoSpect.c statistics.h falsealarm.h TwoSpectSpecFunc.h \
                           IHS.h candidates.h templates.h antenna.h vectormath.h cdfwchisq.h \
                           upperlimits.h cdfdist.h SFTfunctions.h TwoSpectTypes.h TwoSpect.h

lalapps_TwoSpectTemplateBank_SOURCES = TwoSpectTemplateBank.c templates.c vectormath.c \
                                       TwoSpectSpecFunc.c statistics.c templates.h vectormath.h \
                                       TwoSpectSpecFunc.h statistics.h

EXTRA_PROGRAMS = skygridsetup compareCandidates testVectorMath computeSignalDetector

skygridsetup_SOURCES = helperprograms/skygridsetup.c antenna.c antenna.h

compareCandidates_SOURCES = helperprograms/compareCandidates.c

testVectorMath_SOURCES = helperprograms/testVectorMath.c vectormath.c vectormath.h

computeSignalDetector_SOURCES = helperprograms/computeSignalDetector.c TwoSpectSpecFunc.c \
				TwoSpectSpecFunc.h
