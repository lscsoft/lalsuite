AM_CPPFLAGS = -I$(top_builddir)/include -DPKG_DATA_DIR='"$(pkgdatadir)/"'

lalmetaioincludedir = $(includedir)/lal
lalmetaioinclude_HEADERS = \
	LALMetaIOConfig.h \
	LALMetaIOLibVCSInfo.h \
	LALMetaIOVCSInfo.h \
	LIGOLwXML.h \
	LIGOLwXMLArray.h \
	LIGOLwXMLHeaders.h \
	LIGOLwXMLInspiralHeaders.h \
	LIGOLwXMLRead.h \
	LIGOMetadataTables.h \
	LIGOMetadataUtils.h \
	lalmetaioswig_test.i

lib_LTLIBRARIES = liblalmetaio.la
liblalmetaio_la_SOURCES = \
	CreateMetaTableDir.c \
	LALMetaIOVCSInfo.c \
	LIGOLwXML.c \
	LIGOLwXMLArray.c \
	LIGOLwXMLRead.c \
	LIGOMetadataUtils.c

liblalmetaio_la_LDFLAGS = -version-info $(LIBVERSION)

EXTRA_DIST = \
	LALMetaIOConfig.h.in \
	generate_user_env.awk

MOSTLYCLEANFILES = \
	lalmetaio-user-env.shell \
	lalmetaio-user-env.sh \
	lalmetaio-user-env.csh

sysconf_DATA = \
	lalmetaio-user-env.sh \
	lalmetaio-user-env.csh

BUILT_SOURCES = \
	lalmetaio-user-env.sh \
	lalmetaio-user-env.csh

lalmetaio-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# generic user environment script" >$@
	@echo "set LALMETAIO_PREFIX $(prefix)" >>$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@
	@echo "prepend PYTHONPATH $(pyexecdir)" >>$@
	@echo "prepend OCTAVE_PATH $(octexecdir)" >>$@

lalmetaio-user-env.csh lalmetaio-user-env.sh: lalmetaio-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED=$(SED) -v package=$(PACKAGE_NAME) -v output=$@ $< || { rm -f $@; exit 1; }

if HAVE_GIT_REPO
BUILT_SOURCES += LALMetaIOVCSInfo.h
MOSTLYCLEANFILES += LALMetaIOVCSInfo.h
LALMetaIOVCSInfo.h .PHONY : gen-LALMetaIOVCSInfo.h
gen-LALMetaIOVCSInfo.h :
	$(AM_V_at)$(GENERATE_VCS_INFO) LALMetaIOVCSInfo.h $(srcdir)/LALMetaIOVCSInfo.h.git
endif

@LAL_HEADER_LINK@

swig_headers = $(lalmetaioinclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
