## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -DPKG_DATA_DIR='"$(pkgdatadir)/"' @AM_CPPFLAGS@
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_TwoSpect

lalapps_TwoSpect_SOURCES = cmdline.c statistics.c fastchisqinv.c IHS.c candidates.c templates.c antenna.c \
                           vectormath.c cdfwchisq.c upperlimits.c TwoSpect.c \
			   cmdline.h TwoSpectTypes.h statistics.h fastchisqinv.h IHS.h candidates.h templates.h \
                           antenna.h vectormath.h cdfwchisq.h upperlimits.h TwoSpect.h

EXTRA_PROGRAMS = skygridsetup compareCandidates testVectorMath

skygridsetup_SOURCES = helperprograms/skygridsetup.c helperprograms/cmdline_skygridsetup.c antenna.c \
                       helperprograms/cmdline_skygridsetup.h antenna.h

compareCandidates_SOURCES = helperprograms/compareCandidates.c helperprograms/cmdline_compareCandidates.c \
                            helperprograms/cmdline_compareCandidates.h

testVectorMath_SOURCES = helperprograms/testVectorMath.c vectormath.c templates.c candidates.c statistics.c \
                         fastchisqinv.c cdfwchisq.c vectormath.h templates.h candidates.h statistics.h \
                         fastchisqinv.h cdfwchisq.h

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
