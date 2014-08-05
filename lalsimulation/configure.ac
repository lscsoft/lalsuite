AC_PREREQ([2.63])
AC_INIT([LALSimulation],[1.1.0.1],[lal-discuss@ligo.org])
AC_CONFIG_HEADERS([src/config.h src/LALSimulationConfig.h])
AC_CONFIG_SRCDIR([src/LALSimulation.c])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([ \
  lalsimulation.spec \
  lalsimulation.pc \
  lalsimulation-uninstalled.pc \
  Makefile \
  debian/Makefile \
  include/Makefile \
  include/lal/Makefile \
  src/LALSimulationVCSInfo.c \
  src/LALSimulationVCSInfo.h \
  src/Makefile \
  python/Makefile \
  python/lalsimulation/Makefile \
  python/lalsimulation/__init__.py \
  swig/Makefile \
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

# checks for programs
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_AWK
AC_PROG_SED
LALSUITE_CHECK_PYTHON([2.6])

# check for supported mac os x version
LALSUITE_OSX_VERSION_CHECK

# link tests using libtool
LALSUITE_USE_LIBTOOL

# check for header files
AC_HEADER_STDC

# check for gethostname in unistd.h
AC_MSG_CHECKING([for gethostname prototype in unistd.h])
AC_EGREP_HEADER([gethostname],[unistd.h],[AC_MSG_RESULT([yes])]
AC_DEFINE([HAVE_GETHOSTNAME_PROTOTYPE],[1],
          [Define if gethostname prototype is in unistd.h]),
          [AC_MSG_RESULT]([no]))

# check for required libraries
AC_CHECK_LIB([m],[main],,[AC_MSG_ERROR([could not find the math library])])

# check for OpenMP
AC_OPENMP

# check for gsl
PKG_CHECK_MODULES([GSL],[gsl],[true],[false])
LALSUITE_ADD_FLAGS([C],[${GSL_CFLAGS}],[${GSL_LIBS}])
AC_CHECK_LIB([gslcblas],[main],,[AC_MSG_ERROR([could not find the gsl library])])
AC_CHECK_LIB([gsl],[gsl_version],,[AC_MSG_ERROR([could not find the gsl library])])
LALSUITE_CHECK_GSL_VERSION([1.13])
LALSUITE_ENABLE_FAST_GSL

# check for gsl headers
AC_CHECK_HEADERS([gsl/gsl_errno.h],,[AC_MSG_ERROR([could not find the gsl/gsl_errno.h header])])

# check for lal libraries and headers
LALSUITE_CHECK_LIB([LAL],[6.12.0],[LALVersion],[lal/LALStdio.h])
LALSUITE_CHECK_LIB([LALSupport],[6.12.0],[XLALCacheImport],[lal/LALCache.h])

# configure SWIG wrapping modules
LALSUITE_USE_SWIG

# Disable OpenMP by default.
# FIXME: OpenMP should be on by default, but it breaks condor_compiling lalapps.
AS_IF(
  [test "x$enable_openmp" = "x"],
  [OPENMP_CFLAGS=]
)
# add compiler flags for OpenMP
AS_IF(
  [test "x$OPENMP_CFLAGS" = "x"],
  [
    OPENMP_ENABLE_VAL="DISABLED"
    LALSUITE_ADD_FLAGS([C],[-Wno-unknown-pragmas],[])
  ],
  [
    OPENMP_ENABLE_VAL="ENABLED"
    LALSUITE_ADD_FLAGS([C],[${OPENMP_CFLAGS}],[])
  ]
)

# create Automake conditional so that we know whether we need to build
# OpenMP-specific tests
AM_CONDITIONAL([HAVE_OPENMP],[test "x$OPENMP_ENABLE_VAL" = "xENABLED"])

# add gcc specific flags
LALSUITE_ADD_GCC_FLAGS

# link any executables against LALSimulation library
AC_SUBST([LDADD], ["\$(top_builddir)/src/liblalsimulation.la"])

# add extra preprocessor flags
AM_CPPFLAGS="-DPKG_DATA_DIR='\"\$(pkgdatadir)/\"' ${AM_CPPFLAGS}"
AM_CPPFLAGS="-I\$(top_builddir)/include ${AM_CPPFLAGS}"

# force headers to be installed in ${includedir}/lal, ignoring library name
AC_SUBST([pkgincludedir], ["\${includedir}/lal"])

# force all python modules into the exec dir to prevent mixed C and
# pure-python library code from being split across lib and lib64 on
# RedHat-derived systems.
AC_SUBST([pkgpythondir], ["\${pkgpyexecdir}"])

# output
AC_OUTPUT
AC_MSG_NOTICE([
==================================================
LALSimulation has now been successfully configured:

* Python support is $PYTHON_ENABLE_VAL
* SWIG binding generation is $SWIG_GENERATE_ENABLE_VAL
* SWIG bindings for Octave are $SWIG_BUILD_OCTAVE_ENABLE_VAL
* SWIG bindings for Python are $SWIG_BUILD_PYTHON_ENABLE_VAL
* OpenMP acceleration is $OPENMP_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALSimulation,
and run 'make install' to install LALSimulation.
==================================================
])
