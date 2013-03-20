## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_TwoSpect

lalapps_TwoSpect_SOURCES = cmdline.c statistics.c fastchisqinv.c IHS.c candidates.c templates.c antenna.c \
                           vectormath.c cdfwchisq.c upperlimits.c TwoSpect.c \
			   cmdline.h TwoSpectTypes.h statistics.h fastchisqinv.h IHS.h candidates.h templates.h \
                           antenna.h vectormath.h cdfwchisq.h upperlimits.h TwoSpect.h

EXTRA_PROGRAMS = TwoSpect_templateTest skygridsetup compareCandidates

TwoSpect_templateTest_SOURCES = cmdline.c statistics.c fastchisqinv.c IHS.c candidates.c helperprograms/templates_templateTest.c \
                           antenna.c vectormath.c cdfwchisq.c upperlimits.c helperprograms/TwoSpect_templateTest.c \
                           cmdline.h TwoSpectTypes.h statistics.h fastchisqinv.h IHS.h candidates.h templates.h \
                           antenna.h vectormath.h cdfwchisq.h upperlimits.h TwoSpect.h

skygridsetup_SOURCES = helperprograms/skygridsetup.c helperprograms/cmdline_skygridsetup.c antenna.c \
                       helperprograms/cmdline_skygridsetup.h antenna.h

compareCandidates_SOURCES = helperprograms/compareCandidates.c

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
