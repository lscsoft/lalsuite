## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libdate.la
libdate_la_SOURCES = Julian.c \
                     GPStoUTC.c \
                     DateString.c LMST1.c \
                     SecsToLALDate.c date_value.h TimeDelay.c \
                     TimeDelayFromEarthCenter.c GPStoFloat.c \
                     GPStoINT8.c IncrementGPS.c GPSTimeNow.c \
                     Playground.c SkyAnnulus.c StrToGPS.c \
		     XLALTime.c XLALCivilTime.c XLALSiderealTime.c XLALLeapSeconds.h
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
