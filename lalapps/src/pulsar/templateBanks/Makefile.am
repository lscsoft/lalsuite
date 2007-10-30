## Process this file with automake to produce Makefile.in
# $Id$

AM_CFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS =

EXTRA_PROGRAMS = DirectedSearchTemplates

DirectedSearchTemplates_SOURCES = DirectedSearchTemplates.c FlatLatticeTiling.c DirectedSearch.c
