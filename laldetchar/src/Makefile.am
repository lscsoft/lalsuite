AM_CPPFLAGS = -I$(top_builddir)/include -DPKG_DATA_DIR='"$(pkgdatadir)/"'
LDADD = liblaldetchar.la

laldetcharincludedir = $(includedir)/lal
laldetcharinclude_HEADERS = \
	LALDetChar.h \
	LALDetCharConfig.h \
	LALDetCharGlib.h \
	LALDetCharGlibTypes.h \
	LALDetCharLibVCSInfo.h \
	LALDetCharRange.h \
	LALDetCharHvetoUtils.h \
	LALDetCharHveto.h \
	LALDetCharVCSInfo.h \
	laldetcharswig_test.i

lib_LTLIBRARIES = liblaldetchar.la
liblaldetchar_la_SOURCES = \
	LALDetCharGlibTypes.c \
	LALDetCharInspiralRange.c \
	LALDetCharHvetoUtils.c \
	LALDetCharHveto.c \
	LALDetCharVCSInfo.c

liblaldetchar_la_LDFLAGS = -version-info $(LIBVERSION)

MOSTLYCLEANFILES = \
	laldetchar-user-env.shell \
	laldetchar-user-env.sh \
	laldetchar-user-env.csh

sysconf_DATA = \
	laldetchar-user-env.sh \
	laldetchar-user-env.csh

#bin_PROGRAMS = \
	#laldetchar-hveto

#laldetchar_hveto_SOURCES = hveto.c

BUILT_SOURCES = \
	laldetchar-user-env.sh \
	laldetchar-user-env.csh

laldetchar-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# generic user environment script" >$@
	@echo "set LALDETCHAR_PREFIX $(prefix)" >>$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend PYTHONPATH $(pyexecdir)" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@
	@echo "prepend OCTAVE_PATH $(octexecdir)" >>$@

laldetchar-user-env.csh laldetchar-user-env.sh: laldetchar-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED=$(SED) -v package=$(PACKAGE_NAME) -v output=$@ $< || { rm -f $@; exit 1; }

EXTRA_DIST = \
	LALDetCharConfig.h.in \
	config.h.in \
	generate_user_env.awk

if HAVE_GIT_REPO
BUILT_SOURCES += LALDetCharVCSInfo.h
MOSTLYCLEANFILES += LALDetCharVCSInfo.h
LALDetCharVCSInfo.h .PHONY : gen-LALDetCharVCSInfo.h
gen-LALDetCharVCSInfo.h :
	$(AM_V_at)$(GENERATE_VCS_INFO) LALDetCharVCSInfo.h $(srcdir)/LALDetCharVCSInfo.h.git
endif

DISTCLEANFILES = config.h

@LAL_HEADER_LINK@

swig_headers = $(laldetcharinclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
