AC_PREREQ([2.63])
AC_INIT([lal],[6.11.0.1],[lal-discuss@ligo.org])
AC_CONFIG_HEADERS([include/config.h include/lal/LALConfig.h])
AC_CONFIG_SRCDIR([packages/std/include/LALStdio.h])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([ \
  lal.spec \
  lal.pc \
  lalsupport.pc \
  Makefile \
  debian/Makefile \
  gnuscripts/Makefile \
  include/Makefile \
  include/lal/Makefile \
  lib/Makefile \
  lib/lal/Makefile \
  lib/lalsupport/Makefile \
  lib/lalsupport/doc/Makefile \
  lib/lalsupport/include/Makefile \
  lib/lalsupport/src/Makefile \
  lib/lalsupport/test/Makefile \
  python/Makefile \
  python/__init__.py \
  python/utils/Makefile \
  python/spectrum/Makefile \
  swig/Makefile \
  swig/octave/Makefile \
  swig/python/Makefile \
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

LALSUITE_ARG_VAR
LALSUITE_CHECK_GIT_REPO
LALSUITE_DISTCHECK_CONFIGURE_FLAGS
LALSUITE_WITH_CFLAGS_LIBS

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
AC_SUBST([LIBCURRENT],[7])
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
AC_SUBST([LIBCURRENT_SUPPORT],[4])
AC_SUBST([LIBREVISION_SUPPORT],[0])
AC_SUBST([LIBAGE_SUPPORT],[0])
AC_SUBST([LIBVERSION_SUPPORT],[${LIBCURRENT_SUPPORT}:${LIBREVISION_SUPPORT}:${LIBAGE_SUPPORT}])

LAL_ENABLE_INTELFFT

LALSUITE_WITH_CUDA
LALSUITE_WITH_NVCC_CFLAGS

LALSUITE_ENABLE_BOINC
LALSUITE_ENABLE_GCC_FLAGS

LALSUITE_ENABLE_DEBUG
LAL_ENABLE_MACROS
LAL_ENABLE_PTHREAD_LOCK
LALSUITE_ENABLE_NIGHTLY

# enable SWIG wrapping modules
LALSUITE_ENABLE_SWIG

# get version
lal_version_major=`echo "$VERSION" | cut -d. -f1`
lal_version_minor=`echo "$VERSION" | cut -d. -f2`
lal_version_micro=`echo "$VERSION" | cut -d. -f3`
lal_version_devel=`echo "$VERSION" | cut -d. -f4-`
test -z "$lal_version_micro" && lal_version_micro=0
test -z "$lal_version_devel" && lal_version_devel=0
lal_configure_date=`date +"%Y-%m-%dT%H:%M:%S%z"`
AC_DEFINE_UNQUOTED([LAL_VERSION],["$VERSION"],[LAL Version])
AC_DEFINE_UNQUOTED([LAL_VERSION_MAJOR],[$lal_version_major],[LAL Version Major Number])
AC_DEFINE_UNQUOTED([LAL_VERSION_MINOR],[$lal_version_minor],[LAL Version Minor Number])
AC_DEFINE_UNQUOTED([LAL_VERSION_MICRO],[$lal_version_micro],[LAL Version Micro Number])
AC_DEFINE_UNQUOTED([LAL_VERSION_DEVEL],[$lal_version_devel],[LAL Version Devel Number])
AC_DEFINE_UNQUOTED([LAL_CONFIGURE_ARGS],["$ac_configure_args"],[LAL Configure Arguments])
AC_DEFINE_UNQUOTED([LAL_CONFIGURE_DATE],["$lal_configure_date"],[LAL Configure Date])

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
    LALSUITE_ADD_LIBS(${PTHREAD_LIBS})
    LALSUITE_ADD_CFLAGS(${PTHREAD_CFLAGS})
    CC="$PTHREAD_CC"
  else
    AC_MSG_ERROR([do not know how to compile posix threads])
  fi
fi

# check for programs
LALSUITE_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_AWK
AC_PROG_SED
AM_PATH_PYTHON([2.6],,[AC_MSG_ERROR([Python-2.6, or higher, is required])])

# check for supported mac os x version
LALSUITE_OSX_VERSION_CHECK

# link tests using libtool
LALSUITE_USE_LIBTOOL

# set additional compiler flags
AC_MSG_CHECKING([for additional compiler flags])
case "${host_os}" in
  solaris*) LALSUITE_ADD_CFLAGS(-D_POSIX_C_SOURCE=200112L)
  AC_MSG_RESULT([-D_POSIX_C_SOURCE=200112L]);;
  *) AC_MSG_RESULT([none]);;
esac

# initialise lalsupport libs
LALSUPPORT_LIBADD=""

# system library checks
AC_CHECK_LIB([m],[sin])
AC_CHECK_LIB([z],[compress],[LALSUPPORT_LIBADD="$LALSUPPORT_LIBADD -lz"],
             [AC_MSG_ERROR([could not find the zlib library])])

# check for platform specific libs
case "${host_os}" in
  solaris*) AC_CHECK_LIB([sunmath],[sincosp]);;
esac

# check for system headers files
AC_HEADER_STDC
AC_CHECK_HEADERS([sys/time.h unistd.h getopt.h malloc.h regex.h glob.h])
AC_CHECK_HEADERS([zlib.h],,[AC_MSG_ERROR([could not find zlib.h])])
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

# check for gsl
PKG_CHECK_MODULES([GSL],[gsl],[true],[false])
LALSUITE_ADD_CFLAGS(${GSL_CFLAGS})
LALSUITE_ADD_LIBS(${GSL_LIBS})
AC_CHECK_LIB([gslcblas],[main])
AC_CHECK_LIB([gsl],[gsl_strerror],,[AC_MSG_ERROR([could not find the gsl library])])
LAL_CHECK_GSL_VERSION([1.9])
LALSUITE_ENABLE_FAST_GSL

# check for gsl headers
AC_CHECK_HEADERS([gsl/gsl_errno.h],,[AC_MSG_ERROR([could not find the gsl/gsl_errno.h header])])

# check for fft libraries
if test "${intelfft}" = "false" ; then
  fftw3="true"
  PKG_CHECK_MODULES([FFTW3],[fftw3 fftw3f],[true],[false])
  LALSUITE_ADD_CFLAGS(${FFTW3_CFLAGS})
  LALSUITE_ADD_LIBS(${FFTW3_LIBS})
  AC_CHECK_LIB([fftw3f],[fftwf_execute_dft],,[AC_MSG_ERROR([could not find the fftw3f library])],[-lm])
  AC_CHECK_LIB([fftw3],[fftw_execute_dft],,[AC_MSG_ERROR([could not find the fftw3 library])],[-lm])
else
  AC_MSG_WARN([Using Intel FFT routines])
  if test "x${qthread}" = "xtrue" ; then
    LAL_INTEL_MKL_QTHREAD_WARNING
  else
    AC_CHECK_LIB([pthread],[pthread_create],,[AC_MSG_ERROR([could not find the pthread library])])
  fi
  AC_CHECK_LIB([guide],[__kmp_abort],,[AC_MSG_ERROR([could not find the guide library included with Intel MKL])])
  AC_CHECK_HEADERS([pthread.h],,[AC_MSG_ERROR([could not find the pthread.h header])])
  if test "${enable_static}" = "yes"; then
    if test "${enable_shared}" = "yes" ; then
      LAL_INTEL_FFT_LIBS_MSG_ERROR
    else
      AC_CHECK_LIB([mkl_ia32],[DftiCreateDescriptor],,[AC_MSG_ERROR([could not find the Intel FFT library])])
    fi
  else
    AC_CHECK_LIB([mkl],[DftiCreateDescriptor],,[AC_MSG_ERROR([could not find the Intel FFT library])])
  fi
  AC_CHECK_HEADERS([mkl_dfti.h],,[AC_MSG_ERROR([could not find the Intel FFT header])])
fi

# check for fft headers
if test "${intelfft}" = "false" ; then
  AC_CHECK_HEADERS([fftw3.h],,[AC_MSG_ERROR([could not find the fftw3.h header])])
fi

# check for boinc
if test "${boinc}" = "true" ; then
  if test -n "${BOINC_PREFIX}" ; then
    # do compilation checks with c++
    AC_LANG_PUSH([C++])

    # set boinc flags
    BOINC_CFLAGS="-I${BOINC_PREFIX}/include"
    BOINC_LIBS="-L${BOINC_PREFIX}/lib -lboinc_api -lboinc"
    LALSUITE_ADD_CFLAGS(${BOINC_CFLAGS})
    LALSUITE_ADD_LIBS(${BOINC_LIBS})

    # check for boinc libraries
    AC_CHECK_LIB([boinc],[boinc_fopen],[true],[boinc=false])
    AC_CHECK_LIB([boinc_api],[boinc_finish],[true],[boinc=false],[-lboinc])

    # check for boinc headers
    AC_CHECK_HEADERS([boinc/boinc_api.h],,[boinc=false])

    # go back to c
    AC_LANG_POP([C++])
  fi
fi
LALSUITE_ENABLE_MODULE([BOINC],[boinc])

# define what is being enabled
if test "${cuda}" = "true"; then
  AC_DEFINE([LAL_CUDA_ENABLED],[1],[Define if using cuda library])
fi
if test "${fftw3}" = "true"; then
  AC_DEFINE([LAL_FFTW3_ENABLED],[1],[Define if using fftw3 library])
fi
if test "${boinc}" = "true"; then
  AC_DEFINE([LAL_BOINC_ENABLED],[1],[Define if using boinc library])
  AC_SUBST(BOINC_CFLAGS)
  AC_SUBST(BOINC_LIBS)
fi

# configure SWIG wrapping modules
LALSUITE_USE_SWIG([XLAL LAL])

# record the libraries needed for linking lal
LAL_LIBADD="$LIBS"
AC_SUBST(LAL_LIBADD)
AC_SUBST(LALSUPPORT_LIBADD)

# checks for typedefs, structures, and compiler characteristics
AC_C_BIGENDIAN
AC_C_CONST
AC_C_STRINGIZE
AC_TYPE_SIZE_T
AC_HEADER_TIME

# checks for library functions
AC_CHECK_FUNC([gmtime_r],[:],[
  AC_CHECK_FUNCS([gmtime_s],[:],[AC_MSG_ERROR([could not find either gmtime_r() or gmtime_s() functions])])
])
AC_CHECK_FUNC([localtime_r],[:],[
  AC_CHECK_FUNCS([localtime_s],[:],[AC_MSG_ERROR([could not find either localtime_r() or localtime_s() functions])])
])
AC_CHECK_FUNCS([stat],,[AC_MSG_ERROR([could not find the stat() function])])
AC_CHECK_FUNCS([getopt],,[AC_MSG_ERROR([could not find the getopt() function])])
AC_CHECK_FUNCS([getopt_long],,[AC_MSG_ERROR([could not find the getopt_long() function])])
AC_CHECK_FUNCS([putenv])

# set package documentation install path
pkgdocdir="\${prefix}/share/doc/${PACKAGE}-${VERSION}"
AC_SUBST([pkgdocdir])

# set automake conditionals
AM_CONDITIONAL([INTELFFT],[test x$intelfft = xtrue])
AM_CONDITIONAL([QTHREAD],[test x$qthread = xtrue])

# add gcc specific flags --- must be done at end or tests will break!
LALSUITE_ADD_GCC_FLAGS

# force all python modules into the exec dir to prevent mixed C and
# pure-python library code from being split across lib and lib64 on
# RedHat-derived systems.
AC_SUBST([pkgpythondir], ["\${pkgpyexecdir}"])

# ensure CPPFLAGS/CFLAGS and LDFLAGS/LIBS are in correct order
LALSUITE_ADD_CFLAGS()
LALSUITE_ADD_LIBS()

INCLUDEMKFILE=$srcdir/gnuscripts/include.mk
AC_SUBST_FILE([INCLUDEMKFILE])
AC_OUTPUT
echo "
==================================================
LAL has now been successfully configured:

* CUDA support is $CUDA_ENABLE_VAL
* BOINC support is $BOINC_ENABLE_VAL
* SWIG wrapping module for Octave is $SWIG_OCTAVE_ENABLE_VAL
* SWIG wrapping module for Python is $SWIG_PYTHON_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LAL,
and run 'make install' to install LAL.
==================================================
"
