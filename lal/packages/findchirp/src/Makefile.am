## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = liblal.la

if FRAMEANDMPI
liblal_la_SOURCES = FindChirpExch.c
else
liblal_la_SOURCES = 
endif

EXTRA_DIST       = FindChirpExch.c

CC=@MPICC@

INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include
