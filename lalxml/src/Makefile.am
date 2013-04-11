AM_CPPFLAGS = -I$(top_builddir)/include

lalxmlincludedir = $(includedir)/lal
lalxmlinclude_HEADERS = \
	LALXML.h \
	LALXMLConfig.h \
	LALXMLLibVCSInfo.h \
	LALXMLVCSInfo.h \
	LALXMLVOTableCommon.h \
	LALXMLVOTableSerializers.h \
	lalxmlswig_test.i

lib_LTLIBRARIES = liblalxml.la
liblalxml_la_SOURCES = \
	LALXML.c \
	LALXMLVCSInfo.c \
	LALXMLVOTableCommon.c \
	LALXMLVOTableSerializers.c
liblalxml_la_LDFLAGS = -version-info $(LIBVERSION)

MOSTLYCLEANFILES = \
	lalxml-user-env.shell \
	lalxml-user-env.sh \
	lalxml-user-env.csh

sysconf_DATA = \
	lalxml-user-env.sh \
	lalxml-user-env.csh

BUILT_SOURCES = \
	lalxml-user-env.sh \
	lalxml-user-env.csh

lalxml-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# generic user environment script" >$@
	@echo "set LALXML_PREFIX $(prefix)" >>$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend PYTHONPATH $(pyexecdir)" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@
	@printf "$(SWIG_USER_ENV)" >>$@

lalxml-user-env.csh lalxml-user-env.sh: lalxml-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED=$(SED) -v package=$(PACKAGE_NAME) -v output=$@ $< || { rm -f $@; exit 1; }

EXTRA_DIST = \
	generate_user_env.awk \
	LALXMLConfig.h.in

if HAVE_GIT_REPO
BUILT_SOURCES += LALXMLVCSInfo.h
MOSTLYCLEANFILES += LALXMLVCSInfo.h
LALXMLVCSInfo.h .PHONY : gen-LALXMLVCSInfo.h
gen-LALXMLVCSInfo.h :
	$(AM_V_at)$(GENERATE_VCS_INFO) LALXMLVCSInfo.h $(srcdir)/LALXMLVCSInfo.h.git
endif

@LAL_HEADER_LINK@

swig_headers = $(lalxmlinclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
