## Process this file with automake to produce Makefile.in
# $Id$
bin_PROGRAMS = lalapps_findchirp lalapps_inspinj lalapps_findchirp_post \
               lalapps_inspawgfile lalapps_exc_resp
lalapps_findchirp_SOURCES = findchirp.c findchirp.h findchirp_args.c
lalapps_inspinj_SOURCES = inspinj.c
lalapps_findchirp_post_SOURCES = findchirp_post.c
lalapps_inspawgfile_SOURCES = inspawgfile.c
lalapps_exc_resp_SOURCES = exc_resp.c
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
MOSTLYCLEANFILES = injepochs.ilwd injparams.ilwd
