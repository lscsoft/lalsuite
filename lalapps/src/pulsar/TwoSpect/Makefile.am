## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_TwoSpect

lalapps_TwoSpect_SOURCES = cmdline.c IHS.c candidates.c templates.c antenna.c TwoSpect.c \
			   cmdline.h TwoSpectTypes.h IHS.h candidates.h templates.h antenna.h TwoSpect.h
