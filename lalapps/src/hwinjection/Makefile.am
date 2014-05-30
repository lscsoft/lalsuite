LDADD =
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/lalapps.am

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
