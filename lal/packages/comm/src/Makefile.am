## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libcomm.la

if MPI
libcomm_la_SOURCES = Comm.c
else
libcomm_la_SOURCES = 
endif

EXTRA_DIST       = Comm.c

CC=@MPICC@

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

