## Process this file with automake to produce Makefile.in
# $Id$
bin_PROGRAMS = lalapps_findchirp lalapps_inspinj
lalapps_findchirp_SOURCES = findchirp.c findchirp.h findchirp_args.c
lalapps_inspinj_SOURCES = inspinj.c
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
MOSTLYCLEANFILES = injepochs.ilwd injparams.ilwd
