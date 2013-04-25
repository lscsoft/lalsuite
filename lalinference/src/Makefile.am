AM_CPPFLAGS = -I$(top_builddir)/include -DPKG_DATA_DIR='"$(pkgdatadir)/"'

if LALXML
XMLHDR = LALInferenceXML.h
XMLSRC = LALInferenceXML.c
endif

lalinferenceincludedir = $(includedir)/lal
lalinferenceinclude_HEADERS = \
	LALInference.h \
	LALInferenceConfig.h \
	LALInferenceLikelihood.h \
	LALInferenceNestedSampler.h \
	LALInferencePrior.h \
	LALInferenceReadData.h \
	LALInferenceTemplate.h \
	LALInferenceProposal.h \
	LALInferenceLibVCSInfo.h \
	LALInferenceInit.h \
	LALInferenceVCSInfo.h \
        LALInferenceRemoveLines.h \
	lalinferenceswig_test.i \
	$(XMLHDR)

lib_LTLIBRARIES = liblalinference.la
liblalinference_la_SOURCES = \
	LALInference.c \
	LALInferenceLikelihood.c \
	LALInferenceNestedSampler.c \
	LALInferencePrior.c \
	LALInferenceReadData.c \
	LALInferenceTemplate.c \
	LALInferenceProposal.c \
	LALInferenceInitCBC.c \
	LALInferenceVCSInfo.c \
        LALInferenceRemoveLines.c \
	$(XMLSRC)

liblalinference_la_LDFLAGS = -version-info $(LIBVERSION)

MOSTLYCLEANFILES = \
	lalinference-user-env.shell \
	lalinference-user-env.sh \
	lalinference-user-env.csh

sysconf_DATA = \
	lalinference-user-env.sh \
	lalinference-user-env.csh

BUILT_SOURCES = \
	lalinference-user-env.sh \
	lalinference-user-env.csh

lalinference-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# generic user environment script" >$@
	@echo "set LALINFERENCE_PREFIX $(prefix)" >>$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend PYTHONPATH $(pyexecdir)" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@
	@printf "$(SWIG_USER_ENV)" >>$@

lalinference-user-env.csh lalinference-user-env.sh: lalinference-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED=$(SED) -v package=$(PACKAGE_NAME) -v output=$@ $< || { rm -f $@; exit 1; }

EXTRA_DIST = \
	LALInferenceConfig.h.in \
	LALInferenceXML.c \
	LALInferenceXML.h \
	chapter_inference.dox \
	generate_user_env.awk

if HAVE_GIT_REPO
BUILT_SOURCES += LALInferenceVCSInfo.h
MOSTLYCLEANFILES += LALInferenceVCSInfo.h
LALInferenceVCSInfo.h .PHONY : gen-LALInferenceVCSInfo.h
gen-LALInferenceVCSInfo.h :
	$(AM_V_at)$(GENERATE_VCS_INFO) LALInferenceVCSInfo.h $(srcdir)/LALInferenceVCSInfo.h.git
endif

@LAL_HEADER_LINK@

swig_headers = $(lalinferenceinclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
