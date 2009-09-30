INCLUDES = -I$(top_builddir)/include

lalstochasticincludedir = $(includedir)/lal
lalstochasticinclude_HEADERS = \
	AstroOmega.h \
	SimulatePopcorn.h \
	SimulateSB.h \
	StochasticCrossCorrelation.h \
	StochasticMC.h

lib_LTLIBRARIES = liblalstochastic.la
liblalstochastic_la_SOURCES = \
	AstroOmega.c \
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

CLEANFILES = lalstochastic-user-env.sh lalstochastic-user-env.csh
sysconf_DATA = lalstochastic-user-env.sh lalstochastic-user-env.csh
BUILT_SOURCE = lalstochastic-user-env.sh lalstochastic-user-env.csh

lalstochastic-user-env.sh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LAL Stochastic" >$@ ;\
	echo "LALSTOCHASTIC_PREFIX=\"$(prefix)\"" >>$@ ;\
	echo "export LALSTOCHASTIC_PREFIX" >>$@ ;\
	echo "PATH=\"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "MANPATH=\"$(mandir):\$${MANPATH}\"" >>$@ ;\
	echo "LD_LIBRARY_PATH=\"$(libdir):\$${LD_LIBRARY_PATH}\"" >>$@ ;\
	echo "DYLD_LIBRARY_PATH=\"$(libdir):\$${DYLD_LIBRARY_PATH}\"" >>$@ ;\
	echo "PKG_CONFIG_PATH=\"$(libdir)/pkgconfig:\$${PKG_CONFIG_PATH}\"" >>$@ ;\
	echo "export PATH MANPATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH PKG_CONFIG_PATH" >>$@

lalstochastic-user-env.csh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LAL Stochastic" >$@ ;\
	echo "setenv LALSTOCHASTIC_PREFIX \"$(prefix)\"" >>$@ ;\
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

@LAL_HEADER_LINK@
