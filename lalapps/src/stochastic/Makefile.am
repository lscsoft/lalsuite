LDADD =
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/lalapps.am

## sources
STOCH_SRCS       = data_input.c data_input.h \
                   data_output.c data_output.h \
                   sgwb.c sgwb.h misc.c misc.h \
                   stochastic.c stochastic.h
STOCH_FAKE_SRCS  = fake_data.c fake_data.h \
                   data_output.c data_output.h \
                   sgwb.c sgwb.h misc.c misc.h \
                   stochastic_fake.c stochastic.h

STOCH_VIRGO_SRCS = stochastic_virgo.c stochastic_virgo.h
STOCH_PREPROCESS_SRCS = stochastic_preprocess.c stochastic_virgo.h

STOPP_SRCS       = stopp.c
STOPP_BAYES_SRCS = stopp_bayes.c
OLAPREDFCN_SRCS  = olapredfcn.c olapredfcn.h \
                   olapredfcn_args.c

## frame and stochastic library dependant programs
if FRAME
if LALFRAME
if LALSTOCHASTIC
FRAMEPROGRAMS = lalapps_stochastic lalapps_stochastic_fake lalapps_stochastic_virgo lalapps_stochastic_preprocess
lalapps_stochastic_SOURCES      = $(STOCH_SRCS)
lalapps_stochastic_fake_SOURCES = $(STOCH_FAKE_SRCS)
lalapps_stochastic_virgo_SOURCES = $(STOCH_VIRGO_SRCS)
lalapps_stochastic_preprocess_SOURCES = $(STOCH_PREPROCESS_SRCS)
endif
endif
endif

## metaio and stochastic library dependant programs
if LALMETAIO
if LALSTOCHASTIC
METAIOPROGRAMS = lalapps_stopp lalapps_stopp_bayes
lalapps_stopp_SOURCES           = $(STOPP_SRCS)
lalapps_stopp_bayes_SOURCES     = $(STOPP_BAYES_SRCS)
endif
endif

## stochastic library programs
if LALSTOCHASTIC
STOCHASTICPROGRAMS = lalapps_olapredfcn
lalapps_olapredfcn_SOURCES      = $(OLAPREDFCN_SRCS)
endif

bin_PROGRAMS      = ${STOCHASTICPROGRAMS} ${FRAMEPROGRAMS} ${METAIOPROGRAMS}

if LALSTOCHASTIC
if HAVE_PYTHON
bin_SCRIPTS = lalapps_stochastic_pipe lalapps_stochastic_bayes
pkgpython_PYTHON = stochastic.py
endif
endif

if LALSTOCHASTIC
EXTRA_DIST        = stochastic.py stochastic_pipe.in stochastic_bayes.in
endif
MOSTLYCLEANFILES  = *STOCHASTIC*.xml *CCSPECTRA*.gwf
CLEANFILES        = python_config.sed lalapps_stochastic_pipe \
                    lalapps_stochastic_bayes

## rules to build python scripts

python_config.sed:
	@-rm -f python_config.sed
	@echo "s+@PYTHONLIBDIR@+@pkgpythondir@+g" > python_config.sed
	@echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed

lalapps_stochastic_pipe: $(srcdir)/stochastic_pipe.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/stochastic_pipe.in > \
		lalapps_stochastic_pipe
	@chmod +x lalapps_stochastic_pipe

lalapps_stochastic_bayes: $(srcdir)/stochastic_bayes.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/stochastic_bayes.in > \
		lalapps_stochastic_bayes
	@chmod +x lalapps_stochastic_bayes
