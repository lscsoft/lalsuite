## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = liblal.la

if MPI
liblal_la_SOURCES = Comm.c
else
liblal_la_SOURCES = 
endif

EXTRA_DIST       = Comm.c

CC=@MPICC@

INCLUDES = -I$(top_srcdir)/include -I$(srcdir)/../include

