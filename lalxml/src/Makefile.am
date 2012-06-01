INCLUDES = -I$(top_builddir)/include

lalxmlincludedir = $(includedir)/lal
lalxmlinclude_HEADERS = \
	LALXML.h \
	LALXMLConfig.h \
	LALXMLVCSInfo.h \
	LALXMLVOTableCommon.h \
	LALXMLVOTableSerializers.h

lib_LTLIBRARIES = liblalxml.la
liblalxml_la_SOURCES = \
	LALXML.c \
	LALXMLVCSInfo.c \
	LALXMLVOTableCommon.c \
	LALXMLVOTableSerializers.c
liblalxml_la_LDFLAGS = -version-info $(LIBVERSION)

CLEANFILES = \
	lalxml-user-env.shell \
	lalxml-user-env.sh \
	lalxml-user-env.csh

sysconf_DATA = \
	lalxml-user-env.sh \
	lalxml-user-env.csh

BUILT_SOURCES = \
	._vcsinfo \
	LALXMLVCSInfo.h \
	lalxml-user-env.sh \
	lalxml-user-env.csh

lalxml-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# generic user environment script" >$@
	@echo "set LALXML_PREFIX $(prefix)" >>$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend LD_LIBRARY_PATH $(libdir)" >>$@
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@
	@echo "$(SWIG_USER_ENV)" | $(SED) 's|~E~O~L~|\n|g' >>$@

lalxml-user-env.csh lalxml-user-env.sh: lalxml-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED=$(SED) -v package=$(PACKAGE_NAME) -v output=$@ $< || { rm -f $@; exit 1; }

EXTRA_DIST = \
	generate_user_env.awk \
	generate_vcs_info.py \
	LALXMLConfig.h.in \
	LALXMLVCSInfo.h.in

LALXMLVCSInfo.h: ._vcsinfo

._vcsinfo: $(srcdir)/generate_vcs_info.py LALXMLConfig.h
	@$(PYTHON) $(srcdir)/generate_vcs_info.py LALXML $(srcdir) $$PWD

@LAL_HEADER_LINK@

swig_headers = $(lalxmlinclude_HEADERS)
include $(top_srcdir)/swig/swig-header.am
