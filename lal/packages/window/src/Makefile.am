## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libwindow.la
libwindow_la_SOURCES = Window.c

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
