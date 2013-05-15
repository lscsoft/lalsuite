AM_CPPFLAGS = -I$(top_builddir)/include -DPKG_DATA_DIR='"$(pkgdatadir)/"'

lalstochasticincludedir = $(includedir)/lal
lalstochasticinclude_HEADERS = \
	AstroOmega.h \
	LALStochasticConfig.h \
	LALStochasticLibVCSInfo.h \
	LALStochasticVCSInfo.h \
	LIGOLwXMLStochasticRead.h \
	SimulatePopcorn.h \
	SimulateSB.h \
	StochasticCrossCorrelation.h \
	StochasticMC.h \
	lalstochasticswig_test.i

lib_LTLIBRARIES = liblalstochastic.la
liblalstochastic_la_SOURCES = \
	AstroOmega.c \
	LALStochasticVCSInfo.c \
	LIGOLwXMLStochasticRead.c \
	OverlapReductionFunction.c \
	SimulatePopcorn.c \
	SimulateSB.c \
	StochasticCrossCorrelation.c \
	StochasticInverseNoise.c \
	StochasticMC.c \
	StochasticOmegaGW.c \
	StochasticOptimalFilter.c \
	StochasticOptimalFilterNormalization.c \
	ZeroPadAndFFT.c
liblalstochastic_la_LDFLAGS = -version-info $(LIBVERSION)

MOSTLYCLEANFILES = \
	lalstochastic-user-env.shell \
	lalstochastic-user-env.sh \
	lalstochastic-user-env.csh

sysconf_DATA = \
	lalstochastic-user-env.sh \
	lalstochastic-user-env.csh

BUILT_SOURCES = \
	lalstochastic-user-env.sh \
	lalstochastic-user-env.csh

lalstochastic-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# generic user environment script" >$@
	@echo "set LALSTOCHASTIC_PREFIX $(prefix)" >>$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend PYTHONPATH $(pyexecdir)" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@
	@echo "prepend OCTAVE_PATH $(octexecdir)" >>$@

lalstochastic-user-env.csh lalstochastic-user-env.sh: lalstochastic-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED=$(SED) -v package=$(PACKAGE_NAME) -v output=$@ $< || { rm -f $@; exit 1; }

EXTRA_DIST = \
	LALStochasticConfig.h.in \
	config.h.in \
	generate_user_env.awk

DISTCLEANFILES = config.h

if HAVE_GIT_REPO
BUILT_SOURCES += LALStochasticVCSInfo.h
MOSTLYCLEANFILES += LALStochasticVCSInfo.h
LALStochasticVCSInfo.h .PHONY : gen-LALStochasticVCSInfo.h
gen-LALStochasticVCSInfo.h :
	$(AM_V_at)$(GENERATE_VCS_INFO) LALStochasticVCSInfo.h $(srcdir)/LALStochasticVCSInfo.h.git
endif

@LAL_HEADER_LINK@

swig_headers = $(lalstochasticinclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
