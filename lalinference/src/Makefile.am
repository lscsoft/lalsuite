INCLUDES = -I$(top_builddir)/include

if LALXML
XMLHDR = LALInferenceXML.h
XMLSRC = LALInferenceXML.c
endif

lalinferenceincludedir = $(includedir)/lal

HDR= \
	LALInference.h \
	LALInferenceConfig.h \
	LALInferenceLikelihood.h \
	LALInferenceNestedSampler.h \
	LALInferencePrior.h \
	LALInferenceReadData.h \
	LALInferenceTemplate.h \
	LALInferenceVCSInfo.h \
	$(XMLHDR)


lib_LTLIBRARIES = liblalinference.la

SRC= \
	LALInference.c \
	LALInferenceLikelihood.c \
	LALInferenceNestedSampler.c \
	LALInferencePrior.c \
	LALInferenceReadData.c \
	LALInferenceTemplate.c \
	LALInferenceVCSInfo.c \
	$(XMLSRC)

lalinferenceinclude_HEADERS = $(HDR)

liblalinference_la_SOURCES = $(SRC)

liblalinference_la_LDFLAGS = -version-info $(LIBVERSION)

CLEANFILES = \
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

lalinference-user-env.sh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LAL Inference" >$@ ;\
	echo "LALINFERENCE_PREFIX=\"$(prefix)\"" >>$@ ;\
	echo "export LALINFERENCE_PREFIX" >>$@ ;\
	echo "PATH=\"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "MANPATH=\"$(mandir):\$${MANPATH}\"" >>$@ ;\
	echo "LD_LIBRARY_PATH=\"$(libdir):\$${LD_LIBRARY_PATH}\"" >>$@ ;\
	echo "DYLD_LIBRARY_PATH=\"$(libdir):\$${DYLD_LIBRARY_PATH}\"" >>$@ ;\
	echo "PKG_CONFIG_PATH=\"$(libdir)/pkgconfig:\$${PKG_CONFIG_PATH}\"" >>$@ ;\
	echo "export PATH MANPATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH PKG_CONFIG_PATH" >>$@

lalinference-user-env.csh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LAL Inference" >$@ ;\
	echo "setenv LALINFERENCE_PREFIX \"$(prefix)\"" >>$@ ;\
	echo "setenv PATH \"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "if ( \$$?MANPATH ) then" >>$@ ;\
	echo "	setenv MANPATH \"$(mandir):\$${MANPATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "	setenv MANPATH \"$(mandir)\"" >>$@ ;\
	echo "endif" >>$@ ;\
	echo "if ( \$$?LD_LIBRARY_PATH ) then" >>$@ ;\
	echo "	setenv LD_LIBRARY_PATH \"$(libdir):\$${LD_LIBRARY_PATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "	setenv LD_LIBRARY_PATH \"$(libdir)\"" >>$@ ;\
	echo "endif" >>$@ ;\
	echo "if ( \$$?DYLD_LIBRARY_PATH ) then" >>$@ ;\
	echo "	setenv DYLD_LIBRARY_PATH \"$(libdir):\$${DYLD_LIBRARY_PATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "	setenv DYLD_LIBRARY_PATH \"$(libdir)\"" >>$@ ;\
	echo "endif" >>$@ ;\
	echo "if ( \$$?PKG_CONFIG_PATH ) then" >>$@ ;\
	echo "	setenv PKG_CONFIG_PATH \"$(libdir)/pkgconfig:\$${PKG_CONFIG_PATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "	setenv PKG_CONFIG_PATH \"$(libdir)/pkgconfig\"" >>$@ ;\
	echo "endif" >>$@

EXTRA_DIST = \
	LALInferenceConfig.h.in \
	LALInferenceVCSInfo.h.in \
	LALInferenceXML.c \
	LALInferenceXML.h \
	generate_vcs_info.py

LALInferenceVCSInfo.h: ._vcsinfo

._vcsinfo: $(srcdir)/generate_vcs_info.py LALInferenceConfig.h
	@echo "creating LALInferenceVCSInfo.h"
	@$(PYTHON) $(srcdir)/generate_vcs_info.py --sed LALInference $(srcdir) $$PWD

@LAL_HEADER_LINK@
