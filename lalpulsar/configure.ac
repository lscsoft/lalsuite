AC_PREREQ([2.63])
AC_INIT([lalpulsar],[1.9.0.1],[lal-discuss@ligo.org])
AC_CONFIG_HEADERS([src/config.h src/LALPulsarConfig.h])
AC_CONFIG_SRCDIR([src/ComputeFstat.c])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([ \
  lalpulsar.spec \
  lalpulsar.pc \
  lalpulsar-uninstalled.pc \
  Makefile \
  gnuscripts/Makefile \
  debian/Makefile \
  include/Makefile \
  include/lal/Makefile \
  src/Makefile \
  python/Makefile \
  python/__init__.py \
  swig/Makefile \
  swig/octave/Makefile \
  swig/python/Makefile \
  test/Makefile
])
AM_INIT_AUTOMAKE([1.11 foreign subdir-objects color-tests parallel-tests])

AC_CANONICAL_HOST

LALSUITE_CHECK_GIT_REPO
LALSUITE_DISTCHECK_CONFIGURE_FLAGS

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
AC_SUBST([LIBCURRENT],[9])
AC_SUBST([LIBREVISION],[0])
AC_SUBST([LIBAGE],[0])
AC_SUBST([LIBVERSION],[${LIBCURRENT}:${LIBREVISION}:${LIBAGE}])

# nightly build
LALSUITE_ENABLE_NIGHTLY

# enable LALXML support
LALSUITE_ENABLE_LALXML

# boinc
LALSUITE_ENABLE_BOINC

# enable strict gcc flags
LALSUITE_ENABLE_GCC_FLAGS

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

# check for SIMD extensions
LALPULSAR_WITH_SIMD

# checks for programs
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_AWK
AC_PROG_SED
AM_PATH_PYTHON([2.6],,[AC_MSG_ERROR([Python-2.6, or higher, is required])])

# check for supported mac os x version
LALSUITE_OSX_VERSION_CHECK

# link tests using libtool
LALSUITE_USE_LIBTOOL

# check for header files
AC_HEADER_STDC

# check for required libraries
AC_CHECK_LIB([m],[main],,[AC_MSG_ERROR([could not find the math library])])

# check for gsl
PKG_CHECK_MODULES([GSL],[gsl],[true],[false])
LALSUITE_ADD_FLAGS([C],[${GSL_CFLAGS}],[${GSL_LIBS}])
AC_CHECK_LIB([gslcblas],[main],,[AC_MSG_ERROR([could not find the gsl library])])
AC_CHECK_LIB([gsl],[gsl_version],,[AC_MSG_ERROR([could not find the gsl library])])
LALSUITE_CHECK_GSL_VERSION([1.13])
LALSUITE_ENABLE_FAST_GSL

# check for gsl headers
AC_CHECK_HEADERS([gsl/gsl_errno.h],,[AC_MSG_ERROR([could not find the gsl/gsl_errno.h header])])

if test "$lalxml" = "true"; then
  # check for libxml2
  PKG_CHECK_MODULES([LIBXML2],[libxml-2.0 >= 2.6],[true],[false])
  LALSUITE_ADD_FLAGS([C],[${LIBXML2_CFLAGS}],[${LIBXML2_LIBS}])
  AC_CHECK_LIB([xml2],[xmlInitParser],,[AC_MSG_WARN([cannot find the xml2 library])])
fi

# for boinc enabled build do checks for lal in c++ mode
if test "${boinc}" = "true" ; then
  AC_LANG_PUSH([C++])
fi

# check for lal libraries and headers
LALSUITE_CHECK_LIB([LAL],[6.12.0],[LALVersion],[lal/LALStdio.h])

# check that lal was built with BOINC support, if required
LALSUITE_CHECK_BOINC

# check for lalsupport, needed for test suite
LALSUITE_CHECK_LIB([LALSupport],[6.12.0],[XLALCacheImport],[lal/LALCache.h])

# check for lalxml
LALSUITE_CHECK_OPT_LIB([LALXML],[1.1.5],[XLALXMLFilePrintElements],[lal/LALXML.h])

# back to c mode
if test "${boinc}" = "true" ; then
  AC_LANG_POP([C++])
fi

# configure SWIG wrapping modules
LALSUITE_USE_SWIG([XLAL LAL])

# add gcc specific flags
LALSUITE_ADD_GCC_FLAGS

# force all python modules into the exec dir to prevent mixed C and
# pure-python library code from being split across lib and lib64 on
# RedHat-derived systems.
AC_SUBST([pkgpythondir], ["\${pkgpyexecdir}"])

LAL_HEADER_LINK=${ac_aux_dir}/header-link.mk
AC_SUBST_FILE([LAL_HEADER_LINK])

# output
AC_OUTPUT
AC_MSG_NOTICE([
==================================================
LALPulsar has now been successfully configured:

* SIMD extensions are $SIMD_ENABLE_VAL
* LALXML library support is $LALXML_ENABLE_VAL
* BOINC support is $BOINC_ENABLE_VAL
* SWIG wrapping module for Octave is $SWIG_OCTAVE_ENABLE_VAL
* SWIG wrapping module for Python is $SWIG_PYTHON_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALPulsar,
and run 'make install' to install LALPulsar.
==================================================
])
