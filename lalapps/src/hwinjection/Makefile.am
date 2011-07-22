## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if FRAME
frameprogs = lalapps_psinject
else
frameprogs =
endif

bin_PROGRAMS = $(frameprogs)

EXTRA_PROGRAMS = online_psinject

lalapps_psinject_SOURCES = psinject.c

## ---------- compile ONLINE psinject code for LIGO/GEO hardware injections ----------
online_psinject_SOURCES = psinject.c SIStr.c SIStr.h
online_psinject_CPPFLAGS = -DONLINE -I$(GDS_DIR)/include -I$(GDS_DIR)/include/dtt
online_psinject_LDFLAGS = -L$(GDS_DIR)/lib -lawg -ltestpoint
SIStr.c: check_GDS

.PHONY: check_GDS
check_GDS:
	@if test -z "${GDS_DIR}"; then \
	echo "================================================================================";\
	echo "WARNING: Environment-variable 'GDS_DIR' was not set for online_psinject !";\
	echo "To compile this, you can use the syntax";\
	echo "$ GDS_DIR=/your/GDS/dir make online_psinject";\
	echo "================================================================================";\
	exit 1;\
	fi
## --------------------------------------------------------------------------------

TESTS_ENVIRONMENT = \
LALSUITE_TOP_SRCDIR=$(LALSUITE_TOP_SRCDIR) \
LAL_DATA_PATH=.$${LALPULSAR_PREFIX:+":$${LALPULSAR_PREFIX}/share/lalpulsar"}$${LALSUITE_TOP_SRCDIR:+":$${LALSUITE_TOP_SRCDIR}/lalpulsar/test"}

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/misc/vcsID.common
