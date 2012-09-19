AM_CPPFLAGS = -I$(top_builddir)/include

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
	LALInferenceVCSInfo.h \
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
	LALInferenceVCSInfo.c \
	$(XMLSRC)

liblalinference_la_LDFLAGS = -version-info $(LIBVERSION)

CLEANFILES = \
	lalinference-user-env.shell \
	lalinference-user-env.sh \
	lalinference-user-env.csh

sysconf_DATA = \
	lalinference-user-env.sh \
	lalinference-user-env.csh

BUILT_SOURCES = \
	._vcsinfo \
	LALInferenceVCSInfo.h \
	lalinference-user-env.sh \
	lalinference-user-env.csh

lalinference-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# generic user environment script" >$@
	@echo "set LALINFERENCE_PREFIX $(prefix)" >>$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend LD_LIBRARY_PATH $(libdir)" >>$@
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@

lalinference-user-env.csh lalinference-user-env.sh: lalinference-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED=$(SED) -v package=$(PACKAGE_NAME) -v output=$@ $< || { rm -f $@; exit 1; }

EXTRA_DIST = \
	LALInferenceConfig.h.in \
	LALInferenceVCSInfo.h.in \
	LALInferenceXML.c \
	LALInferenceXML.h \
	generate_user_env.awk \
	generate_vcs_info.py

LALInferenceVCSInfo.h: ._vcsinfo

._vcsinfo: $(srcdir)/generate_vcs_info.py LALInferenceConfig.h
	@$(PYTHON) $(srcdir)/generate_vcs_info.py LALInference $(srcdir) $$PWD

@LAL_HEADER_LINK@
