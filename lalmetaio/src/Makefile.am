AM_CPPFLAGS = -I$(top_builddir)/include

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
	LIGOMetadataUtils.h

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
	LALMetaIOVCSInfo.h.in \
	generate_user_env.awk \
	generate_vcs_info.py

CLEANFILES = \
	lalmetaio-user-env.shell \
	lalmetaio-user-env.sh \
	lalmetaio-user-env.csh

sysconf_DATA = \
	lalmetaio-user-env.sh \
	lalmetaio-user-env.csh

BUILT_SOURCES = \
	._vcsinfo \
	LALMetaIOVCSInfo.h \
	lalmetaio-user-env.sh \
	lalmetaio-user-env.csh

lalmetaio-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# generic user environment script" >$@
	@echo "set LALMETAIO_PREFIX $(prefix)" >>$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend LD_LIBRARY_PATH $(libdir)" >>$@
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@

lalmetaio-user-env.csh lalmetaio-user-env.sh: lalmetaio-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED=$(SED) -v package=$(PACKAGE_NAME) -v output=$@ $< || { rm -f $@; exit 1; }

LALMetaIOVCSInfo.h: ._vcsinfo

._vcsinfo: $(srcdir)/generate_vcs_info.py LALMetaIOConfig.h
	@$(PYTHON) $(srcdir)/generate_vcs_info.py LALMetaIO $(srcdir) $$PWD

@LAL_HEADER_LINK@
