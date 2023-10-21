AC_PREREQ([2.63])
AC_INIT([LALPulsar],[6.0.6.1],[lal-discuss@ligo.org])
AC_CONFIG_HEADERS([lib/config.h lib/LALPulsarConfig.h])
AC_CONFIG_SRCDIR([lib/ComputeFstat.c])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])

# set library api info
# 1. start with version information of "0:0:0" for each libtool library
# 2. update the version information only immediately before a public
#    release of your software. More frequent updates are unnecessary,
#    and only guarantee that the current interface number gets larger
#    faster.
# 3. if the library source code has changed at all since the last
#    update, then increment revision ("c:r:a" becomes "c:r+1:a").
# 4. if any interfaces have been added, removed, or changed since the
#    last update, increment current, and set revision to 0.
# 5. if any interfaces have been added since the last public release,
#    then increment age.
# 6. if any interfaces have been removed since the last public release,
#    then set age to 0.
AC_SUBST([LIBCURRENT],[29])
AC_SUBST([LIBREVISION],[3])
AC_SUBST([LIBAGE],[0])
AC_SUBST([LIBVERSION],[${LIBCURRENT}:${LIBREVISION}:${LIBAGE}])
AS_CASE([${host_os}],
  [freebsd-aout|qnx|sco|sunos], [AC_SUBST([LIBMAJOR],[${LIBCURRENT}])],
  [AC_SUBST([LIBMAJOR],[$((LIBCURRENT-LIBAGE))])]
)

AC_CONFIG_FILES([ \
  Makefile \
  conda/Makefile \
  conda/meta.yaml.in \
  conda/migrations/Makefile \
  debian/Makefile \
  debian/control \
  debian/liblalpulsar${LIBMAJOR}.install:debian/liblalpulsar.install.in \
  doxygen/Makefile \
  include/Makefile \
  include/lal/Makefile \
  lalpulsar.pc \
  lalpulsar.spec \
  python/Makefile \
  python/lalpulsar/Makefile \
  python/lalpulsar/__init__.py \
  lib/LALPulsarVCSInfo.c \
  lib/LALPulsarVCSInfo.h \
  lib/Makefile \
  bin/CreateEphemeris/Makefile \
  bin/CrossCorr/Makefile \
  bin/FITSTools/Makefile \
  bin/Fscan/Makefile \
  bin/Fstatistic/Makefile \
  bin/GCT/Makefile \
  bin/HWInjection/Makefile \
  bin/HeterodyneSearch/Makefile \
  bin/Hough/Makefile \
  bin/HoughFstat/Makefile \
  bin/longtests/Makefile \
  bin/MakeData/Makefile \
  bin/Makefile \
  bin/SFTTools/Makefile \
  bin/SidebandSearch/Makefile \
  bin/Tools/Makefile \
  bin/TwoSpect/Makefile \
  bin/Weave/Makefile \
  bin/version.c \
  swig/Makefile \
  test/Makefile \
  test/python/Makefile \
  test/python/test_python_imports.py \
  test/python/test_python_scripts.py \
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
MIN_LAL_VERSION="7.4.0"
MIN_LALFRAME_VERSION="3.0.0"
MIN_LALSIMULATION_VERSION="5.3.0"
MIN_LALINFERENCE_VERSION="4.1.0"
AC_SUBST([MIN_LAL_VERSION])
AC_SUBST([MIN_LALFRAME_VERSION])
AC_SUBST([MIN_LALSIMULATION_VERSION])
AC_SUBST([MIN_LALINFERENCE_VERSION])

AC_CANONICAL_HOST

LALSUITE_ARG_VAR

LALSUITE_CHECK_PLATFORM
LALSUITE_CHECK_GIT_REPO
LALSUITE_DISTCHECK_CONFIGURE_FLAGS
LALSUITE_WITH_FALLBACK_DATA_PATH

# nightly build
LALSUITE_ENABLE_NIGHTLY

# enable cfitsio support
LALSUITE_ENABLE_CFITSIO

# enable SIStr support
LALPULSAR_ENABLE_SISTR

# enable strict gcc flags
LALSUITE_ENABLE_GCC_FLAGS

# enable strict definitions
LALSUITE_ENABLE_STRICT_DEFS

# optional lal libraries
LALSUITE_ENABLE_LALFRAME
LALSUITE_ENABLE_LALSIMULATION
LALSUITE_ENABLE_LALINFERENCE

# enable SWIG wrapping modules
LALSUITE_ENABLE_SWIG

# enable Doxygen documentation
LALSUITE_ENABLE_DOXYGEN

# enable CUDA
LALSUITE_WITH_CUDA
if test "${cuda}" = "true"; then
  AC_DEFINE([LALPULSAR_CUDA_ENABLED],[1],[Define if using cuda library])
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
LALPULSAR_CHECK_ALTIVEC

# choose which ephemeris files to install
LALPULSAR_WITH_EPHEM

# checks for programs
AC_PROG_INSTALL
AC_PROG_MKDIR_P
AC_PROG_LN_S
AC_PROG_AWK
AC_PROG_SED
LALSUITE_CHECK_PYTHON([3.5])

# check for supported mac os x version
LALSUITE_OSX_VERSION_CHECK

# link tests using libtool
LALSUITE_USE_LIBTOOL

# check for header files
AC_CHECK_HEADERS([unistd.h])

# check for specific functions
AC_FUNC_STRNLEN

# check for required libraries
AC_CHECK_LIB([m],[main],,[AC_MSG_ERROR([could not find the math library])])

# check for OpenMP
LALSUITE_ENABLE_OPENMP

# check for gsl
PKG_CHECK_MODULES([GSL],[gsl],[true],[false])
LALSUITE_ADD_FLAGS([C],[${GSL_CFLAGS}],[${GSL_LIBS}])
AC_CHECK_LIB([gsl],[gsl_version],,[AC_MSG_ERROR([could not find the gsl library])])
LALSUITE_CHECK_GSL_VERSION([1.13])
LALSUITE_ENABLE_FAST_GSL

# check for gsl headers
AC_CHECK_HEADERS([gsl/gsl_errno.h],,[AC_MSG_ERROR([could not find the gsl/gsl_errno.h header])])

# check for fftw
PKG_CHECK_MODULES([FFTW3],[fftw3 fftw3f],[true],[false])
LALSUITE_ADD_FLAGS([C],[${FFTW3_CFLAGS}],[${FFTW3_LIBS}])
AC_CHECK_LIB([fftw3f],[fftwf_execute_dft],,[AC_MSG_ERROR([could not find the fftw3f library])],[-lm])
AC_CHECK_LIB([fftw3],[fftw_execute_dft],,[AC_MSG_ERROR([could not find the fftw3 library])],[-lm])

# check for fft headers
AC_CHECK_HEADERS([fftw3.h],,[AC_MSG_ERROR([could not find the fftw3.h header])])

# check for cfitsio
LALSUITE_USE_CFITSIO

# check for SIStr library, for pulsar hardware injections
if test "${sistr}" = "true"; then
  AC_CHECK_HEADERS([gds/dtt/SIStr.h],[break],[AC_MSG_ERROR([could not find the SIStr header])])
  AC_CHECK_LIB([SIStr],[SIStrAppInfo],,[AC_MSG_ERROR([could not find the SIStr library])])
  AC_CHECK_LIB([SIStr],[SIStrErrorMsg],,[AC_MSG_ERROR([could not find the SIStr library])])
  AC_CHECK_LIB([SIStr],[SIStrOpen],,[AC_MSG_ERROR([could not find the SIStr library])])
  AC_CHECK_LIB([SIStr],[SIStrAppend],,[AC_MSG_ERROR([could not find the SIStr library])])
  AC_CHECK_LIB([SIStr],[SIStrClose],,[AC_MSG_ERROR([could not find the SIStr library])])
  AC_SUBST([SISTR_LIBS],["-lSIStr"])
fi
LALSUITE_ENABLE_MODULE([SISTR])

# check for lal libraries and headers
LALSUITE_CHECK_LIB([LAL],[${MIN_LAL_VERSION}],[XLALPerror],[lal/XLALError.h])
LALSUITE_CHECK_LIB([LALSupport],[${MIN_LAL_VERSION}],[XLALCacheImport],[lal/LALCache.h])
LALSUITE_CHECK_OPT_LIB([LALFrame],[${MIN_LALFRAME_VERSION}],[XLALFrStreamOpen],[lal/LALFrStream.h])
LALSUITE_CHECK_OPT_LIB([LALSimulation],[${MIN_LALSIMULATION_VERSION}],[XLALSimDetectorStrainREAL8TimeSeries],[lal/LALSimulation.h])
LALSUITE_CHECK_OPT_LIB([LALInference],[${MIN_LALINFERENCE_VERSION}],[LALInferenceAddVariable],[lal/LALInference.h])

# configure SWIG wrapping modules
LALSUITE_USE_SWIG

# configure Doxygen documentation
LALSUITE_USE_DOXYGEN

# configure CUDA
LALSUITE_USE_CUDA

# add gcc specific flags
LALSUITE_ADD_GCC_FLAGS

# link any executables against LALPulsar library
AC_SUBST([LDADD], ["\$(top_builddir)/lib/liblalpulsar.la"])

# add extra preprocessor flags
AM_CPPFLAGS="-I\$(top_builddir)/include ${AM_CPPFLAGS}"

# force headers to be installed in ${includedir}/lal, ignoring library name
AC_SUBST([pkgincludedir], ["\${includedir}/lal"])

# force all python modules into the exec dir to prevent mixed C and
# pure-python library code from being split across lib and lib64 on
# RedHat-derived systems.
AC_SUBST([pkgpythondir], ["\${pkgpyexecdir}"])

# make sure local library module is in python path
LAL_PYTHON_PATH="\$(abs_top_builddir)/python:${LAL_PYTHON_PATH}"

# output
AC_OUTPUT
AC_MSG_NOTICE([
==================================================
LALPulsar has now been successfully configured:

* Python support is $PYTHON_ENABLE_VAL
* LALFrame library support is $LALFRAME_ENABLE_VAL
* LALSimulation library support is $LALSIMULATION_ENABLE_VAL
* LALInference library support is $LALINFERENCE_ENABLE_VAL
* CFITSIO library support is $CFITSIO_ENABLE_VAL
* SIStr support is $SISTR_ENABLE_VAL
* CUDA support is $CUDA_ENABLE_VAL
* SWIG bindings for Octave are $SWIG_BUILD_OCTAVE_ENABLE_VAL
* SWIG bindings for Python are $SWIG_BUILD_PYTHON_ENABLE_VAL
* OpenMP acceleration is $OPENMP_ENABLE_VAL
* Doxygen documentation is $DOXYGEN_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALPulsar,
and run 'make install' to install LALPulsar.
==================================================
])
