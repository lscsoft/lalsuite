BUILT_SOURCES =
sysconf_DATA =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_library_top.am

ACLOCAL_AMFLAGS = -I gnuscripts

EXTRA_DIST += \
	lalapps.spec \
	$(DATA)

SUBDIRS = \
	debian \
	src \
	man \
	doxygen \
	$(END_OF_LIST)

user_environment = \
	set $(PACKAGE)_PREFIX $(prefix); \
	set $(PACKAGE)_DATADIR $(pkgdatadir); \
	prepend MANPATH $(mandir); \
	prepend PATH $(bindir); \
	prepend PYTHONPATH $(pyexecdir) $(pythondir); \
	$(END_OF_LIST)
