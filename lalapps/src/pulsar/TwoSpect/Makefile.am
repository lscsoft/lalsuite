## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_TwoSpect

lalapps_TwoSpect_SOURCES = cmdline.c statistics.c IHS.c candidates.c templates.c antenna.c cdfwchisq.c TwoSpect.c \
			   cmdline.h TwoSpectTypes.h statistics.h IHS.h candidates.h templates.h antenna.h \
                           cdfwchisq.h TwoSpect.h

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/misc/vcsID.common
