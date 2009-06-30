DATAFIND_SRCS = online_datafind.c

if FRAME
if LALFRAME
FRAME_PROGS = lalapps_online_datafind
lalapps_online_datafind_SOURCES = online_datafind.c
endif
endif

bin_PROGRAMS = ${FRAME_PROGS}

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
