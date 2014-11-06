AC_PREREQ([2.63])
AC_INIT([LAL],[6.13.2.1],[lal-discuss@ligo.org])
AC_CONFIG_HEADERS([include/config.h lib/lal/LALConfig.h])
AC_CONFIG_SRCDIR([packages/std/include/LALStdio.h])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([ \
  lal.spec \
  lal.pc \
  lalsupport.pc \
  Makefile \
  debian/Makefile \
  include/Makefile \
  include/lal/Makefile \
  lib/Makefile \
  lib/lal/LALVCSInfo.c \
  lib/lal/LALVCSInfo.h \
  lib/lal/Makefile \
  lib/lalsupport/Makefile \
  lib/lalsupport/doc/Makefile \
  lib/lalsupport/include/Makefile \
  lib/lalsupport/src/Makefile \
  lib/lalsupport/test/Makefile \
  python/Makefile \
  python/lal/Makefile \
  python/lal/__init__.py \
  python/lal/utils/Makefile \
  python/lal/spectrum/Makefile \
  swig/Makefile \
  packages/Makefile \
  packages/std/Makefile \
  packages/std/doc/Makefile \
  packages/std/include/Makefile
  packages/std/src/Makefile \
  packages/std/test/Makefile \
  packages/tools/Makefile \
  packages/tools/doc/Makefile \
  packages/tools/include/Makefile \
  packages/tools/src/Makefile \
  packages/tools/test/Makefile \
  packages/factories/Makefile \
  packages/factories/doc/Makefile \
  packages/factories/include/Makefile \
  packages/factories/src/Makefile \
  packages/factories/test/Makefile \
  packages/vectorops/Makefile \
  packages/vectorops/doc/Makefile \
  packages/vectorops/include/Makefile \
  packages/vectorops/src/Makefile \
  packages/vectorops/test/Makefile \
  packages/utilities/Makefile \
  packages/utilities/doc/Makefile \
  packages/utilities/include/Makefile \
  packages/utilities/src/Makefile \
  packages/utilities/test/Makefile \
  packages/stats/Makefile \
  packages/stats/doc/Makefile \
  packages/stats/include/Makefile \
  packages/stats/src/Makefile \
  packages/stats/test/Makefile \
  packages/inject/Makefile \
  packages/inject/doc/Makefile \
  packages/inject/include/Makefile \
  packages/inject/src/Makefile \
  packages/inject/test/Makefile \
  packages/date/Makefile \
  packages/date/doc/Makefile \
  packages/date/include/Makefile \
  packages/date/src/Makefile \
  packages/date/test/Makefile \
  packages/tdfilter/Makefile \
  packages/tdfilter/doc/Makefile \
  packages/tdfilter/include/Makefile \
  packages/tdfilter/src/Makefile \
  packages/tdfilter/test/Makefile \
  packages/window/Makefile \
  packages/window/doc/Makefile \
  packages/window/include/Makefile \
  packages/window/src/Makefile \
  packages/window/test/Makefile \
  packages/fft/Makefile \
  packages/fft/doc/Makefile \
  packages/fft/include/Makefile \
  packages/fft/src/Makefile \
  packages/fft/test/Makefile \
  packages/noisemodels/Makefile \
  packages/noisemodels/doc/Makefile \
  packages/noisemodels/include/Makefile \
  packages/noisemodels/src/Makefile \
  packages/noisemodels/test/Makefile
])
AM_INIT_AUTOMAKE([1.11 foreign subdir-objects color-tests parallel-tests])

AC_CANONICAL_HOST

LALSUITE_CHECK_GIT_REPO
LALSUITE_DISTCHECK_CONFIGURE_FLAGS

# set library api info - liblal
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
AC_SUBST([LIBCURRENT],[9])
AC_SUBST([LIBREVISION],[0])
AC_SUBST([LIBAGE],[0])
AC_SUBST([LIBVERSION],[${LIBCURRENT}:${LIBREVISION}:${LIBAGE}])

# set library api info - liblalsupport
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
AC_SUBST([LIBCURRENT_SUPPORT],[5])
AC_SUBST([LIBREVISION_SUPPORT],[0])
AC_SUBST([LIBAGE_SUPPORT],[1])
AC_SUBST([LIBVERSION_SUPPORT],[${LIBCURRENT_SUPPORT}:${LIBREVISION_SUPPORT}:${LIBAGE_SUPPORT}])

LAL_ENABLE_INTELFFT
LAL_ENABLE_FFTW3_MEMALIGN

LALSUITE_WITH_CUDA
LALSUITE_WITH_NVCC_CFLAGS

LALSUITE_ENABLE_BOINC
LALSUITE_ENABLE_GCC_FLAGS

LAL_ENABLE_DEBUG
LAL_ENABLE_MACROS
LAL_ENABLE_PTHREAD_LOCK
LALSUITE_ENABLE_NIGHTLY

# enable SWIG wrapping modules
LALSUITE_ENABLE_SWIG

# define version/configure info
LALSUITE_VERSION_CONFIGURE_INFO

# workaround broken libtool multilib support on sl6
LALSUITE_MULTILIB_LIBTOOL_HACK

# use silent build rules
AM_SILENT_RULES([yes])

# check for required compilers
LALSUITE_PROG_COMPILERS

# check for pthreads
lal_pthread_available="yes"
AX_PTHREAD(,[lal_pthread_available="no"])
if test "$lal_pthread_lock" = "true" ; then
  if test "$lal_pthread_available" = "yes" ; then
    LALSUITE_ADD_FLAGS([C],[${PTHREAD_CFLAGS}],[${PTHREAD_LIBS}])
    CC="$PTHREAD_CC"
  else
    AC_MSG_ERROR([do not know how to compile posix threads])
  fi
fi

# check for programs
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_AWK
AC_PROG_SED
LALSUITE_CHECK_PYTHON([2.6])

# check for supported mac os x version
LALSUITE_OSX_VERSION_CHECK

# link tests using libtool
LALSUITE_USE_LIBTOOL

# set additional compiler flags
AC_MSG_CHECKING([for additional compiler flags])
case "${host_os}" in
  solaris*) LALSUITE_ADD_FLAGS([C],[-D_POSIX_C_SOURCE=200112L],[])
  AC_MSG_RESULT([-D_POSIX_C_SOURCE=200112L]);;
  *) AC_MSG_RESULT([none]);;
esac

# system library checks
AC_CHECK_LIB([m],[sin])

# check for platform specific libs
case "${host_os}" in
  solaris*) AC_CHECK_LIB([sunmath],[sincosp]);;
esac

# check for system headers files
AC_HEADER_STDC
AC_CHECK_HEADERS([sys/time.h sys/resource.h unistd.h getopt.h malloc.h regex.h glob.h])
AC_CHECK_HEADERS([stdint.h],,[AC_MSG_ERROR([could not find stdint.h])])
AC_CHECK_HEADERS([inttypes.h],,[AC_MSG_ERROR([could not find inttypes.h])])
AC_MSG_CHECKING([for gethostname prototype in unistd.h])
AC_EGREP_HEADER([gethostname],[unistd.h],AC_MSG_RESULT([yes])
AC_DEFINE([HAVE_GETHOSTNAME_PROTOTYPE],[1],[Define if gethostname prototype is in unistd.h]),AC_MSG_RESULT([no]))

# check for pthread headers
if test "${lal_pthread_lock}" = "true"; then
  AC_CHECK_HEADERS([pthread.h],[break])
fi

# check for platform specific headers
case "${host_os}" in
  solaris*) AC_CHECK_HEADERS([sunmath.h]);;
esac

# check for zlib libraries and headers
PKG_CHECK_MODULES([ZLIB],[zlib],[true],[false])
LALSUITE_PUSH_UVARS
LALSUITE_ADD_FLAGS([C],[${ZLIB_CFLAGS}],[${ZLIB_LIBS}])
AC_SEARCH_LIBS([compress],[],[:],[AC_MSG_ERROR([could not find the zlib library])])
AC_CHECK_HEADER([zlib.h],[:],[AC_MSG_ERROR([could not find the zlib.h header])])
LALSUITE_POP_UVARS

# check for gsl
PKG_CHECK_MODULES([GSL],[gsl],[true],[false])
LALSUITE_ADD_FLAGS([C],[${GSL_CFLAGS}],[${GSL_LIBS}])
AC_CHECK_LIB([gslcblas],[main],,[AC_MSG_ERROR([could not find the gsl library])])
AC_CHECK_LIB([gsl],[gsl_version],,[AC_MSG_ERROR([could not find the gsl library])])
LALSUITE_CHECK_GSL_VERSION([1.13])
LALSUITE_ENABLE_FAST_GSL

# check for gsl headers
AC_CHECK_HEADERS([gsl/gsl_errno.h],,[AC_MSG_ERROR([could not find the gsl/gsl_errno.h header])])

# check for fft libraries
if test "${intelfft}" = "false" ; then
  fftw3="true"
  PKG_CHECK_MODULES([FFTW3],[fftw3 fftw3f],[true],[false])
  LALSUITE_ADD_FLAGS([C],[${FFTW3_CFLAGS}],[${FFTW3_LIBS}])
  AC_CHECK_LIB([fftw3f],[fftwf_execute_dft],,[AC_MSG_ERROR([could not find the fftw3f library])],[-lm])
  AC_CHECK_LIB([fftw3],[fftw_execute_dft],,[AC_MSG_ERROR([could not find the fftw3 library])],[-lm])
else
  AC_MSG_WARN([Using Intel FFT routines])
  if test "x${qthread}" = "xtrue" ; then
    LAL_INTEL_MKL_QTHREAD_WARNING
  else
    AC_CHECK_LIB([pthread],[pthread_create],,[AC_MSG_ERROR([could not find the pthread library])])
  fi
  # AC_CHECK_LIB([guide],[__kmp_abort],,[AC_MSG_ERROR([could not find the guide library included with Intel MKL])])
  AC_CHECK_HEADERS([pthread.h],,[AC_MSG_ERROR([could not find the pthread.h header])])
  if test "${enable_static}" = "yes"; then
    if test "${enable_shared}" = "yes" ; then
      LAL_INTEL_FFT_LIBS_MSG_ERROR
    else
      AC_CHECK_LIB([mkl_rt],[DftiCreateDescriptor],,[AC_MSG_ERROR([could not find the Intel FFT library])])
    fi
  else
    AC_CHECK_LIB([mkl_rt],[DftiCreateDescriptor],,[AC_MSG_ERROR([could not find the Intel FFT library])])
  fi
  AC_CHECK_HEADERS([mkl_dfti.h],,[AC_MSG_ERROR([could not find the Intel FFT header])])
fi

# check for fft headers
if test "${intelfft}" = "false" ; then
  AC_CHECK_HEADERS([fftw3.h],,[AC_MSG_ERROR([could not find the fftw3.h header])])
fi

# check for boinc
if test "x${boinc}" = xtrue; then
  # do compilation checks with c++
  AC_LANG_PUSH([C++])

  # set boinc flags
  AC_ARG_VAR([BOINC_PREFIX],[BOINC installation directory (optional)])
  AC_ARG_VAR([BOINC_EXTRA_CFLAGS],[extra BOINC compiler flags (optional)])
  AC_ARG_VAR([BOINC_EXTRA_LIBS],[extra BOINC libraries (optional)])
  AC_SUBST(BOINC_CFLAGS,"${BOINC_EXTRA_CFLAGS}")
  AC_SUBST(BOINC_LIBS,"-lboinc_api -lboinc ${BOINC_EXTRA_LIBS}")
  if test "x${BOINC_PREFIX}" != x; then
    BOINC_CFLAGS="-I${BOINC_PREFIX}/include ${BOINC_CFLAGS}"
    BOINC_LIBS="-L${BOINC_PREFIX}/lib ${BOINC_LIBS}"
  fi
  LALSUITE_ADD_FLAGS([C],[${BOINC_CFLAGS}],[${BOINC_LIBS}])

  # check for boinc libraries
  AC_CHECK_LIB([boinc],[boinc_fopen],[:],[AC_MSG_ERROR([could not find the boinc library])])
  AC_CHECK_LIB([boinc_api],[boinc_finish],[:],[AC_MSG_ERROR([could not find the boinc_api library])])

  # check for boinc headers
  AC_CHECK_HEADERS([boinc/boinc_api.h],[:],[AC_MSG_ERROR([could not find the boinc_api.h header])])

  # go back to c
  AC_LANG_POP([C++])
fi

# define what is being enabled
if test "${cuda}" = "true"; then
  AC_DEFINE([LAL_CUDA_ENABLED],[1],[Define if using cuda library])
fi
if test "${fftw3}" = "true"; then
  AC_DEFINE([LAL_FFTW3_ENABLED],[1],[Define if using fftw3 library])
  if test "${fftw3_memalign}" = "true"; then
    AC_DEFINE([LAL_FFTW3_MEMALIGN_ENABLED],[1],[Define if using fftw3 library])
  fi
fi
if test "${boinc}" = "true"; then
  AC_DEFINE([LAL_BOINC_ENABLED],[1],[Define if using boinc library])
  AC_SUBST(BOINC_CFLAGS)
  AC_SUBST(BOINC_LIBS)
fi

# configure SWIG wrapping modules
LALSUITE_USE_SWIG

# checks for typedefs, structures, and compiler characteristics
AC_C_BIGENDIAN
AC_C_CONST
AC_C_STRINGIZE
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_HEADER_TIME

# checks for library functions
AC_CHECK_FUNCS([gmtime_r localtime_r])
AC_CHECK_FUNCS([stat],,[AC_MSG_ERROR([could not find the stat() function])])
AC_CHECK_FUNCS([getopt],,[AC_MSG_ERROR([could not find the getopt() function])])
AC_CHECK_FUNCS([getopt_long],,[AC_MSG_ERROR([could not find the getopt_long() function])])
AC_CHECK_FUNCS([putenv])
AC_CHECK_FUNCS([posix_memalign])

# set package documentation install path
pkgdocdir="\${prefix}/share/doc/${PACKAGE}-${VERSION}"
AC_SUBST([pkgdocdir])

# set automake conditionals
AM_CONDITIONAL([INTELFFT],[test x$intelfft = xtrue])
AM_CONDITIONAL([QTHREAD],[test x$qthread = xtrue])

# add gcc specific flags --- must be done at end or tests will break!
LALSUITE_ADD_GCC_FLAGS

# link any executables against LAL and LALSupport libraries
AC_SUBST([LDADD], ["\$(top_builddir)/lib/lal/liblal.la \$(top_builddir)/lib/lalsupport/src/liblalsupport.la"])

# add extra preprocessor flags
AM_CPPFLAGS="-DPKG_DATA_DIR='\"\$(pkgdatadir)/\"' ${AM_CPPFLAGS}"
AM_CPPFLAGS="-I\$(top_builddir)/include ${AM_CPPFLAGS}"

# force all python modules into the exec dir to prevent mixed C and
# pure-python library code from being split across lib and lib64 on
# RedHat-derived systems.
AC_SUBST([pkgpythondir], ["\${pkgpyexecdir}"])

# output
AC_OUTPUT
AC_MSG_NOTICE([
==================================================
LAL has now been successfully configured:

* Python support is $PYTHON_ENABLE_VAL
* CUDA support is $CUDA_ENABLE_VAL
* BOINC support is $BOINC_ENABLE_VAL
* SWIG binding generation is $SWIG_GENERATE_ENABLE_VAL
* SWIG bindings for Octave are $SWIG_BUILD_OCTAVE_ENABLE_VAL
* SWIG bindings for Python are $SWIG_BUILD_PYTHON_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LAL,
and run 'make install' to install LAL.
==================================================
])
