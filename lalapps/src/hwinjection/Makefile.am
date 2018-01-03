include $(top_srcdir)/gnuscripts/lalapps.am

progs = lalapps_psinject

if GDS
progs += online_psinject
endif

bin_PROGRAMS = $(progs)

lalapps_psinject_SOURCES = psinject.c

## ---------- compile ONLINE psinject code for LIGO/GEO hardware injections ----------
online_psinject_SOURCES = psinject.c SIStr.c SIStr.h
online_psinject_CPPFLAGS = -DONLINE
online_psinject_LDFLAGS = -lawg -ltestpoint
## --------------------------------------------------------------------------------
