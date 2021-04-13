BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalapps_pulsar_test.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

bin_PROGRAMS = lalapps_psinject

if GDS
bin_PROGRAMS += online_psinject
endif

lalapps_psinject_SOURCES = psinject.c

## ---------- compile ONLINE psinject code for LIGO/GEO hardware injections ----------
online_psinject_SOURCES = psinject.c SIStr.c SIStr.h
online_psinject_CPPFLAGS = -DONLINE
online_psinject_LDFLAGS = -lawg -ltestpoint
## --------------------------------------------------------------------------------
