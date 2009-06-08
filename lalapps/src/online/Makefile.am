CACHE_SRCS = online_cache.c

if FRAME
if LALFRAME
FRAME_PROGS = lalapps_online_cache
lalapps_online_cache_SOURCES = online_cache.c
endif
endif

bin_PROGRAMS = ${FRAME_PROGS}

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
