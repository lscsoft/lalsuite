AC_PREREQ([2.63])
AC_INIT([LALApps],[10.1.0.1],[lal-discuss@ligo.org])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_SRCDIR([src/lalapps/version2.c.in])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([\
  Makefile \
  conda/Makefile \
  conda/meta.yaml.in \
  conda/migrations/Makefile \
  debian/Makefile \
  debian/control \
  doxygen/Makefile \
  lalapps.spec \
  man/Makefile \
  src/Makefile \
  src/findchirp/Makefile \
  src/frametools/Makefile \
  src/inspiral/Makefile \
  src/lalapps/__init__.py \
  src/lalapps/LALAppsVCSInfo.c \
  src/lalapps/LALAppsVCSInfo.h \
  src/lalapps/Makefile \
  src/lalapps/version2.c \
  src/power/Makefile \
  src/pulsar/Makefile \
  src/ring/Makefile \
  src/string/Makefile \
])
AC_USE_SYSTEM_EXTENSIONS
AM_INIT_AUTOMAKE([ \
  1.11 \
  foreign \
  subdir-objects \
  color-tests \
  parallel-tests \
  tar-pax \
  dist-xz \
  no-dist-gzip \
])

# version requirements
MIN_FRAMEL_VERSION="8.39.2"
MIN_LAL_VERSION="7.7.0"
MIN_LALFRAME_VERSION="3.0.0"
MIN_LALMETAIO_VERSION="4.0.0"
MIN_LALSIMULATION_VERSION="6.2.0"
MIN_LALBURST_VERSION="2.0.0"
MIN_LALINSPIRAL_VERSION="5.0.0"
MIN_LALINFERENCE_VERSION="4.1.0"
MIN_LALPULSAR_VERSION="7.1.0"
MIN_METAIO_VERSION="8.4.0"
MIN_PYTHON_IGWN_LIGOLW_VERSION="2.1.0"
AC_SUBST([MIN_FRAMEL_VERSION])
AC_SUBST([MIN_LAL_VERSION])
AC_SUBST([MIN_LALFRAME_VERSION])
AC_SUBST([MIN_LALMETAIO_VERSION])
AC_SUBST([MIN_LALSIMULATION_VERSION])
AC_SUBST([MIN_LALBURST_VERSION])
AC_SUBST([MIN_LALINSPIRAL_VERSION])
AC_SUBST([MIN_LALINFERENCE_VERSION])
AC_SUBST([MIN_LALPULSAR_VERSION])
AC_SUBST([MIN_METAIO_VERSION])
AC_SUBST([MIN_PYTHON_IGWN_LIGOLW_VERSION])

LT_INIT()

AC_CANONICAL_HOST

LALSUITE_ARG_VAR

LALSUITE_CHECK_PLATFORM
LALSUITE_CHECK_GIT_REPO
LALSUITE_DISTCHECK_CONFIGURE_FLAGS
LALSUITE_WITH_FALLBACK_DATA_PATH

LALSUITE_ENABLE_GCC_FLAGS
LALSUITE_ENABLE_STRICT_DEFS
LALSUITE_ENABLE_NIGHTLY
LALAPPS_ENABLE_CONDOR
LALAPPS_ENABLE_STATIC_BINARIES

LALAPPS_ENABLE_FFTW
LALAPPS_ENABLE_FRAMEL
LALSUITE_ENABLE_CFITSIO

LALSUITE_ENABLE_LALFRAME
LALSUITE_ENABLE_LALMETAIO
LALSUITE_ENABLE_LALSIMULATION
LALSUITE_ENABLE_LALBURST
LALSUITE_ENABLE_LALINSPIRAL
LALSUITE_ENABLE_LALINFERENCE
LALSUITE_ENABLE_LALPULSAR

# enable Doxygen documentation
LALSUITE_ENABLE_DOXYGEN

# enable CUDA
LALSUITE_WITH_CUDA
if test "${cuda}" = "true"; then
  AC_DEFINE([LALAPPS_CUDA_ENABLED],[1],[Define if using cuda library])
fi

# enable help2man documentation
LALSUITE_ENABLE_HELP2MAN

# define version/configure info
LALSUITE_VERSION_CONFIGURE_INFO

# use silent build rules
AM_SILENT_RULES([yes])

# check for required compilers
LALSUITE_PROG_COMPILERS

# check for SIMD extensions
LALSUITE_CHECK_SIMD

# link tests using libtool
if test "${static_binaries}" = "true"; then
  lalsuite_libtool_flags="-all-static"
fi
LALSUITE_USE_LIBTOOL

# checks for programs
AC_PROG_INSTALL
AC_PROG_MKDIR_P
AC_PROG_LN_S
AC_PROG_AWK
AC_PROG_SED
LALSUITE_CHECK_PYTHON([])   # use default minimum Python version

# check for pager programs
LALSUITE_CHECK_PAGER

# check for condor_compile
if test "x$condor" = "xtrue"; then
  AC_CHECK_PROGS([CONDOR_COMPILE],[condor_compile],[false])
  if test "x$CONDOR_COMPILE" = "xfalse"; then
    AC_MSG_ERROR([condor_compile not found])
  fi

  # does condor_compile suffer from condorbug 19343
  AC_MSG_CHECKING([if condor_compile suffers from bug 19343])
  if $CONDOR_COMPILE $CC -print-prog-name=ld 1>/dev/null 2>/dev/null; then
    AC_MSG_RESULT([no])
  else
    AC_MSG_RESULT([yes])
    LD=/usr/bin/ld
  fi

  CC="$CONDOR_COMPILE $CC"

  CONDOR_ENABLE_VAL="ENABLED"
else
  CONDOR_ENABLE_VAL="DISABLED"
fi

# check for supported mac os x version
LALSUITE_OSX_VERSION_CHECK

# check for system libraries
AC_CHECK_LIB([m],[sin])

# check for system headers
AC_CHECK_HEADERS([unistd.h glob.h])

# check for specific functions
AC_CHECK_FUNC([strdup], [], [AC_MSG_ERROR([could not find the strdup function])])

# check for gsl
PKG_CHECK_MODULES([GSL],[gsl],[true],[false])
LALSUITE_ADD_FLAGS([C],[${GSL_CFLAGS}],[${GSL_LIBS}])
AC_CHECK_LIB([gsl],[gsl_version],,[AC_MSG_ERROR([could not find the gsl library])])
LALSUITE_CHECK_GSL_VERSION([1.13])
LALSUITE_ENABLE_FAST_GSL

# check for gsl headers
AC_CHECK_HEADERS([gsl/gsl_errno.h],,[AC_MSG_ERROR([could not find the gsl/gsl_errno.h header])])

# check for fftw libraries
if test "${fftw}" = "true"; then
  PKG_CHECK_MODULES([FFTW3],[fftw3 fftw3f],[true],[false])
  LALSUITE_ADD_FLAGS([C],[${FFTW3_CFLAGS}],[${FFTW3_LIBS}])
  AC_CHECK_LIB([fftw3f],[fftwf_execute_dft],,[AC_MSG_WARN([could not find the fftw3f library])],[-lm])
  AC_CHECK_LIB([fftw3],[fftw_execute_dft],,[AC_MSG_WARN([could not find the fftw3 library])],[-lm])
fi

# check for fftw headers
if test "${fftw}" = "true"; then
  AC_CHECK_HEADERS([fftw3.h],,[AC_MSG_WARN([could not find the fftw3.h header])])
fi
LALSUITE_ENABLE_MODULE([FFTW])

# check for FrameL
if test "${framel}" = "true"; then
  PKG_CHECK_MODULES([FRAMEL],[framel],[true],[false])
  lal_pre_framel_LIBS="$LIBS"
  LALSUITE_ADD_FLAGS([C],[],[${FRAMEL_LIBS}])
  AC_SEARCH_LIBS([FrLibIni],[framel],,
                 [AC_MSG_WARN([could not find the FrameL library])]
                 [framel="false"]
                 [LIBS="$lal_pre_framel_LIBS"])
fi

# check for FrameL headers
if test "${framel}" = "true"; then
  lal_pre_framel_CPPFLAGS="$CPPFLAGS"
  LALSUITE_ADD_FLAGS([C],[${FRAMEL_CFLAGS}],[])
  AC_CHECK_HEADERS([FrameL.h],,
                   [AC_MSG_WARN([could not find the FrameL.h header])]
                   [framel="false"]
                   [CPPFLAGS="$lal_pre_framel_CPPFLAGS"])
fi
LALSUITE_ENABLE_MODULE([FRAMEL])

# metaio
PKG_CHECK_MODULES([METAIO],[libmetaio >= ${MIN_METAIO_VERSION}],[true],[false])
LALSUITE_ADD_FLAGS([C],[${METAIO_CFLAGS}],[${METAIO_LIBS}])
AC_CHECK_LIB([metaio],[MetaioOpen],,[AC_MSG_ERROR([could not find the metaio library])])
AC_CHECK_HEADERS([metaio.h],,[AC_MSG_ERROR([could not find the metaio.h header])])

# check for cfitsio
LALSUITE_USE_CFITSIO

# check for lal libararies and headers
LALSUITE_CHECK_LIB([LAL],[${MIN_LAL_VERSION}],[XLALPerror],[lal/XLALError.h])
LALSUITE_CHECK_LIB([LALSupport],[${MIN_LAL_VERSION}],[XLALCacheImport],[lal/LALCache.h])
LALSUITE_CHECK_OPT_LIB([LALFrame],[${MIN_LALFRAME_VERSION}],[XLALFrStreamOpen],[lal/LALFrStream.h])
LALSUITE_CHECK_OPT_LIB([LALMetaIO],[${MIN_LALMETAIO_VERSION}],[XLALLIGOLwHasTable],[lal/LIGOLwXMLRead.h])
LALSUITE_CHECK_OPT_LIB([LALSimulation],[${MIN_LALSIMULATION_VERSION}],[XLALSimDetectorStrainREAL8TimeSeries],[lal/LALSimulation.h])
LALSUITE_CHECK_OPT_LIB([LALBurst],[${MIN_LALBURST_VERSION}],[XLALEPSearch],[lal/EPSearch.h])
LALSUITE_CHECK_OPT_LIB([LALInspiral],[${MIN_LALINSPIRAL_VERSION}],[LALInspiralParameterCalc],[lal/LALInspiral.h])
LALSUITE_CHECK_OPT_LIB([LALInference],[${MIN_LALINFERENCE_VERSION}],[LALInferenceAddVariable],[lal/LALInference.h])
LALSUITE_CHECK_OPT_LIB([LALPulsar],[${MIN_LALPULSAR_VERSION}],[XLALComputeFstat],[lal/ComputeFstat.h])

# configure Doxygen documentation
LALSUITE_USE_DOXYGEN

# configure CUDA
LALSUITE_USE_CUDA

# check for gethostname prototype
AC_MSG_CHECKING([for gethostname prototype in unistd.h])
AC_EGREP_HEADER([gethostname],[unistd.h],[AC_MSG_RESULT([yes])]
  [AC_DEFINE([HAVE_GETHOSTNAME_PROTOTYPE],[1],
             [Define if gethostname prototype is in unistd.h])],
  [AC_MSG_RESULT([no])])

AC_CHECK_FUNCS([setenv])

# set prefix
AC_DEFINE_UNQUOTED([PREFIX],["$prefix"],[Install prefix])

LALSUITE_ADD_GCC_FLAGS

# add extra preprocessor flags
AM_CPPFLAGS="-I\$(top_builddir)/src/lalapps ${AM_CPPFLAGS}"
AM_CPPFLAGS="-I\$(top_srcdir)/src/lalapps ${AM_CPPFLAGS}"
AM_CPPFLAGS="-I\$(top_srcdir)/src ${AM_CPPFLAGS}"

# force all python modules into the exec dir to prevent mixed C and
# pure-python library code from being split across lib and lib64 on
# RedHat-derived systems.
AC_SUBST([pkgpythondir], ["\${pkgpyexecdir}"])

# output
AC_OUTPUT
AC_MSG_NOTICE([
==================================================
LALApps has now been successfully configured:

* Python support is $PYTHON_ENABLE_VAL
* FFTW library support is $FFTW_ENABLE_VAL
* FrameL library support is $FRAMEL_ENABLE_VAL
* LALFrame library support is $LALFRAME_ENABLE_VAL
* LALMetaIO library support is $LALMETAIO_ENABLE_VAL
* LALSimulation library support is $LALSIMULATION_ENABLE_VAL
* LALBurst library support is $LALBURST_ENABLE_VAL
* LALInspiral library support is $LALINSPIRAL_ENABLE_VAL
* LALInference library support is $LALINFERENCE_ENABLE_VAL
* LALPulsar legacy support is $LALPULSAR_ENABLE_VAL
* CFITSIO library support is $CFITSIO_ENABLE_VAL
* Condor support is $CONDOR_ENABLE_VAL
* CUDA support is $CUDA_ENABLE_VAL
* Doxygen documentation is $DOXYGEN_ENABLE_VAL
* help2man documentation is $HELP2MAN_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALApps,
and run 'make install' to install LALApps.
==================================================
])
