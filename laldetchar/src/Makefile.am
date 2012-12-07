AM_CPPFLAGS = -I$(top_builddir)/include

laldetcharincludedir = $(includedir)/lal
laldetcharinclude_HEADERS = \
	LALDetChar.h \
	LALDetCharConfig.h \
	LALDetCharLibVCSInfo.h \
	LALDetCharRange.h \
	LALDetCharVCSInfo.h

lib_LTLIBRARIES = liblaldetchar.la
liblaldetchar_la_SOURCES = \
	LALDetCharInspiralRange.c \
	LALDetCharVCSInfo.c

liblaldetchar_la_LDFLAGS = -version-info $(LIBVERSION)

MOSTLYCLEANFILES = \
	laldetchar-user-env.shell \
	laldetchar-user-env.sh \
	laldetchar-user-env.csh

sysconf_DATA = \
	laldetchar-user-env.sh \
	laldetchar-user-env.csh

BUILT_SOURCES = \
	laldetchar-user-env.sh \
	laldetchar-user-env.csh

laldetchar-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# generic user environment script" >$@
	@echo "set LALDETCHAR_PREFIX $(prefix)" >>$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend LD_LIBRARY_PATH $(libdir)" >>$@
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend PYTHONDIR $(pyexecdir)" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@

laldetchar-user-env.csh laldetchar-user-env.sh: laldetchar-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED=$(SED) -v package=$(PACKAGE_NAME) -v output=$@ $< || { rm -f $@; exit 1; }

EXTRA_DIST = \
	LALDetCharConfig.h.in \
	generate_user_env.awk

if HAVE_GIT_REPO
BUILT_SOURCES += LALDetCharVCSInfo.h
MOSTLYCLEANFILES += LALDetCharVCSInfo.h
LALDetCharVCSInfo.h .PHONY : gen-LALDetCharVCSInfo.h
gen-LALDetCharVCSInfo.h :
	$(GENERATE_VCS_INFO) LALDetCharVCSInfo.h $(srcdir)/LALDetCharVCSInfo.h.git
endif

@LAL_HEADER_LINK@

swig_headers = $(laldetcharinclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
