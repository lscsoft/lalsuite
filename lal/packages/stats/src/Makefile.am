include $(top_srcdir)/gnuscripts/lallibs_header_links.am

AM_CPPFLAGS = -I$(top_builddir)/include -I$(srcdir)/../include @AM_CPPFLAGS@

noinst_LTLIBRARIES = libstats.la

libstats_la_SOURCES = \
	LALMoment.c \
	LALCorrelation.c

noinst_HEADERS = LALMoment_source.c
