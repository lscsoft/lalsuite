BUILT_SOURCES =
sysconf_DATA =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_library_top.am

ACLOCAL_AMFLAGS = -I gnuscripts

EXTRA_DIST += \
	lalsimulation.spec \
	$(DATA)

SUBDIRS = \
	debian \
	include \
	src \
	swig \
	python \
	test \
	doxygen \
	$(END_OF_LIST)

user_environment = \
	set $(PACKAGE)_PREFIX $(prefix); \
	set $(PACKAGE)_DATADIR $(pkgdatadir); \
	prepend MANPATH $(mandir); \
	prepend OCTAVE_PATH $(octexecdir); \
	prepend PATH $(bindir); \
	prepend PKG_CONFIG_PATH $(libdir)/pkgconfig; \
	prepend PYTHONPATH $(pyexecdir) $(pythondir); \
	$(END_OF_LIST)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = lalsimulation.pc
