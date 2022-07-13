BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalpulsar_bin.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

if GDS

bin_PROGRAMS = lalpulsar_hwinject

lalpulsar_hwinject_SOURCES = \
	hwinject.c \
	$(END_OF_LIST)

lalpulsar_hwinject_CPPFLAGS = $(AM_CPPFLAGS) -DONLINE
lalpulsar_hwinject_LIBS = $(GDS_LIBS)

else

noinst_PROGRAMS = lalpulsar_hwinject

lalpulsar_hwinject_CPPFLAGS = $(AM_CPPFLAGS)

lalpulsar_hwinject_SOURCES = \
	hwinject.c \
	$(END_OF_LIST)

endif

# Add shell test scripts to this variable
test_scripts += test_hwinject.sh
