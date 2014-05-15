include $(top_srcdir)/gnuscripts/lallibs_header_links.am

AM_CPPFLAGS = -I$(top_builddir)/include -I$(srcdir)/../include @AM_CPPFLAGS@

noinst_LTLIBRARIES = libstd.la
libstd_la_SOURCES = \
	Grid.c \
	LALDebugLevel.c \
	LALError.c \
	LALGSL.c \
	LALMalloc.c \
	LALString.c \
	LALVCSInfoType.c \
	LALVersion.c \
	StringConvert.c \
	StringToken.c \
	XLALError.c

noinst_HEADERS = Grid_source.c
