## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libdate.la
libdate_la_SOURCES = GPStoUTC.c LMST1.c date_value.h TimeDelay.c \
                     GPStoFloat.c IncrementGPS.c GPSTimeNow.c \
                     Playground.c SkyAnnulus.c StrToGPS.c \
		     XLALTime.c XLALCivilTime.c XLALSiderealTime.c XLALLeapSeconds.h
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
