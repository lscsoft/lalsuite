AC_PREREQ([2.63])
AC_INIT([LALFrame],[3.0.7.1],[lal-discuss@ligo.org])
AC_CONFIG_HEADERS([lib/config.h lib/LALFrameConfig.h])
AC_CONFIG_SRCDIR([lib/LALFrameIO.c])
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
AC_SUBST([LIBCURRENT],[14])
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
  debian/liblalframe${LIBMAJOR}.install:debian/liblalframe.install.in \
  doxygen/Makefile \
  include/Makefile \
  include/lal/Makefile \
  lalframe.pc \
  lalframe.spec \
  python/Makefile \
  python/lalframe/Makefile \
  python/lalframe/__init__.py \
  python/lalframe/utils/Makefile \
  lib/LALFrameVCSInfo.c \
  lib/LALFrameVCSInfo.h \
  lib/Makefile \
  bin/Makefile \
  bin/version.c \
  swig/Makefile \
  test/Makefile \
  test/python/Makefile \
  test/python/test_python_imports.py \
  test/python/test_python_scripts.py \
])
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
MIN_LAL_VERSION="7.7.0"
MIN_FRAMECPPC_VERSION="2.5.5"
MIN_FRAMEL_VERSION="8.42.2"
AC_SUBST([MIN_LAL_VERSION])
AC_SUBST([MIN_FRAMECPPC_VERSION])
AC_SUBST([MIN_FRAMEL_VERSION])

AC_CANONICAL_HOST

LALSUITE_ARG_VAR

LALSUITE_CHECK_PLATFORM
LALSUITE_CHECK_GIT_REPO
LALSUITE_DISTCHECK_CONFIGURE_FLAGS
LALSUITE_WITH_FALLBACK_DATA_PATH

# nightly build
LALSUITE_ENABLE_NIGHTLY

# enable check for framec
LALFRAME_ENABLE_FRAMEC

# enable check for framel
LALFRAME_ENABLE_FRAMEL

# enable strict gcc flags
LALSUITE_ENABLE_GCC_FLAGS

# enable strict definitions
LALSUITE_ENABLE_STRICT_DEFS

# enable SWIG wrapping modules
LALSUITE_ENABLE_SWIG

# enable Doxygen documentation
LALSUITE_ENABLE_DOXYGEN

# enable help2man documentation
LALSUITE_ENABLE_HELP2MAN

# define version/configure info
LALSUITE_VERSION_CONFIGURE_INFO

# use silent build rules
AM_SILENT_RULES([yes])

# check for required compilers
LALSUITE_PROG_COMPILERS

# check for pthread, needed for low latency data test codes
AX_PTHREAD([lalframe_pthread=true],[lalframe_pthread=false])
AM_CONDITIONAL([PTHREAD],[test x$lalframe_pthread = xtrue])

# checks for programs
AC_PROG_INSTALL
AC_PROG_MKDIR_P
AC_PROG_LN_S
AC_PROG_AWK
AC_PROG_SED
LALSUITE_CHECK_PYTHON([])   # use default minimum Python version

# check for supported mac os x version
LALSUITE_OSX_VERSION_CHECK

# link tests using libtool
LALSUITE_USE_LIBTOOL

# check for header files
AC_CHECK_HEADERS([unistd.h])

# check for gethostname in unistd.h
AC_MSG_CHECKING([for gethostname prototype in unistd.h])
AC_EGREP_HEADER([gethostname],[unistd.h],[AC_MSG_RESULT([yes])]
AC_DEFINE([HAVE_GETHOSTNAME_PROTOTYPE],[1],
          [Define if gethostname prototype is in unistd.h]),
          [AC_MSG_RESULT]([no]))

# check for required libraries
AC_CHECK_LIB([m],[main],,[AC_MSG_ERROR([could not find the math library])])

# checks for library functions
AC_CHECK_FUNCS([gmtime_r localtime_r])

# check for framec or libframe libraries and headers
PKG_PROG_PKG_CONFIG
FRAMEC_AVAILABLE="no"
FRAMEL_AVAILABLE="no"

AS_IF([test "${framec}" = "true"], [
       FRAMEC_AVAILABLE="yes"
       PKG_CHECK_MODULES([FRAMECPPC],[framecppc >= ${MIN_FRAMECPPC_VERSION}],[true],[false])
       LALSUITE_ADD_FLAGS([C],[${FRAMECPPC_CFLAGS}],[${FRAMECPPC_LIBS}])
       AC_CHECK_LIB([framecppc],[FrameCFileOpen],,[FRAMEC_AVAILABLE="no"])
       AC_CHECK_HEADERS([framecppc/FrameC.h],,[FRAMEC_AVAILABLE="no"])
])

AS_IF([test "${framel}" = "true"], [
       FRAMEL_AVAILABLE="yes"
       PKG_CHECK_MODULES([FRAMEL],[framel >= ${MIN_FRAMEL_VERSION}],[true],[false])
       LALSUITE_ADD_FLAGS([C],[${FRAMEL_CFLAGS}],[${FRAMEL_LIBS}])
       AC_CHECK_LIB([framel],[FrLibIni],,[FRAMEL_AVAILABLE="no"])
       AC_CHECK_HEADERS([FrameL.h],,[FRAMEL_AVAILABLE="no"])
])

# make sure a frame library was found
AS_IF([test "${FRAMEC_AVAILABLE}" = "no" -a "${FRAMEL_AVAILABLE}" = "no"], [
       AC_MSG_ERROR([no frame library available])
])

AS_IF([test "${FRAMEC_AVAILABLE}" = "yes"],
      [AC_SUBST([PC_REQUIRES_FRAMEC], [", framecppc >= ${MIN_FRAMECPPC_VERSION}"])],
      [AC_SUBST([PC_REQUIRES_FRAMEC], [""])]
)

AS_IF([test "${FRAMEL_AVAILABLE}" = "yes"],
      [AC_SUBST([PC_REQUIRES_FRAMEL], [", framel"])],
      [AC_SUBST([PC_REQUIRES_FRAMEL], [""])]
)

AM_CONDITIONAL([FRAMEC],[test "${FRAMEC_AVAILABLE}" = "yes"])
AM_CONDITIONAL([FRAMEL],[test "${FRAMEL_AVAILABLE}" = "yes"])

# check for lal libraries and headers
LALSUITE_CHECK_LIB([LAL],[${MIN_LAL_VERSION}],[XLALPerror],[lal/XLALError.h])
LALSUITE_CHECK_LIB([LALSupport],[${MIN_LAL_VERSION}],[XLALCacheImport],[lal/LALCache.h])

# configure SWIG wrapping modules
LALSUITE_USE_SWIG

# configure Doxygen documentation
LALSUITE_USE_DOXYGEN

# add gcc specific flags
LALSUITE_ADD_GCC_FLAGS

# link any executables against LALFrame library
AC_SUBST([LDADD], ["\$(top_builddir)/lib/liblalframe.la"])

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
LALFrame has now been successfully configured:

* Python support is $PYTHON_ENABLE_VAL
* FrameC availability... ${FRAMEC_AVAILABLE}
* FrameL availability... ${FRAMEL_AVAILABLE}
* SWIG bindings for Octave are $SWIG_BUILD_OCTAVE_ENABLE_VAL
* SWIG bindings for Python are $SWIG_BUILD_PYTHON_ENABLE_VAL
* Doxygen documentation is $DOXYGEN_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALFrame,
and run 'make install' to install LALFrame.
==================================================
])
