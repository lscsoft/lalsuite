## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libfindchirp.la

if FRAMEANDMPI
libfindchirp_la_SOURCES = FindChirpExch.c
else
libfindchirp_la_SOURCES = 
endif

EXTRA_DIST       = FindChirpExch.c

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
