## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if LALFRAME
frameprogs = lalapps_psinject
if GDS
frameprogs += online_psinject
endif
endif

bin_PROGRAMS = $(frameprogs)

lalapps_psinject_SOURCES = psinject.c

## ---------- compile ONLINE psinject code for LIGO/GEO hardware injections ----------
online_psinject_SOURCES = psinject.c SIStr.c SIStr.h
online_psinject_CPPFLAGS = -DONLINE
online_psinject_LDFLAGS = -lawg -ltestpoint
## --------------------------------------------------------------------------------

TESTS_ENVIRONMENT = \
LALSUITE_TOP_SRCDIR=$(LALSUITE_TOP_SRCDIR) \
LAL_DATA_PATH=.$${LALPULSAR_PREFIX:+":$${LALPULSAR_PREFIX}/share/lalpulsar"}$${LALSUITE_TOP_SRCDIR:+":$${LALSUITE_TOP_SRCDIR}/lalpulsar/test"}

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/misc/vcsID.common
