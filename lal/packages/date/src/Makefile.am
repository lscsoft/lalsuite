## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libdate.la
libdate_la_SOURCES = Julian.c UtoGPS.c Utime.c DateString.c LMST1.c \
                     SecsToLALDate.c date_value.h TimeDelay.c \
                     TimeDelayFromEarthCenter.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
