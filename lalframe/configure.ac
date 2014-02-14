AC_PREREQ([2.63])
AC_INIT([lalframe],[1.1.1.1],[lal-discuss@ligo.org])
AC_CONFIG_HEADERS([src/config.h src/LALFrameConfig.h])
AC_CONFIG_SRCDIR([src/LALFrameIO.c])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([ \
  lalframe.spec \
  lalframe.pc \
  lalframe-uninstalled.pc \
  Makefile \
  gnuscripts/Makefile \
  debian/Makefile \
  include/Makefile \
  include/lal/Makefile \
  src/Makefile \
  python/Makefile \
  python/__init__.py \
  python/utils/Makefile \
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
AC_SUBST([LIBCURRENT],[7])
AC_SUBST([LIBREVISION],[1])
AC_SUBST([LIBAGE],[0])
AC_SUBST([LIBVERSION],[${LIBCURRENT}:${LIBREVISION}:${LIBAGE}])

# nightly build
LALSUITE_ENABLE_NIGHTLY

# enable check for inotify
LALFRAME_ENABLE_INOTIFY

# enable strict gcc flags
LALSUITE_ENABLE_GCC_FLAGS

# enable SWIG wrapping modules
LALSUITE_ENABLE_SWIG

# allow use of either FrameC or FrameL
LALFRAME_WITH_FRAME_LIBRARY

# define version/configure info
LALSUITE_VERSION_CONFIGURE_INFO

# workaround broken libtool multilib support on sl6
LALSUITE_MULTILIB_LIBTOOL_HACK

# use silent build rules
AM_SILENT_RULES([yes])

# check for required compilers
LALSUITE_PROG_COMPILERS

# check for pythread, needed for low latency data test codes
AX_PTHREAD([lalframe_pthread=true],[lalframe_pthread=false])
AM_CONDITIONAL([PTHREAD],[test x$lalframe_pthread = xtrue])

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
AS_IF([test "${framec}" = "true"],[
       PKG_CHECK_MODULES([FRAMECPPC],[framecppc >= 1.19.32],[true],[false])
       LALSUITE_ADD_FLAGS([C],[${FRAMECPPC_CFLAGS}],[${FRAMECPPC_LIBS}])
       AC_CHECK_LIB([framecppc],[FrameCFileOpen],,[AC_MSG_ERROR([could not find the framec library])])
       AC_CHECK_HEADERS([framecppc/FrameC.h],,[AC_MSG_ERROR([could not find the FrameC.h header])])
       AC_DEFINE([USE_FRAMEC],[],[Define to use FrameC])
       FRAME_LIBRARY="FrameC"])
AS_IF([test "${framel}" = "true"],[
       PKG_CHECK_MODULES([LIBFRAME],[libframe],[true],[false])
       LALSUITE_ADD_FLAGS([C],[${LIBFRAME_CFLAGS}],[${LIBFRAME_LIBS}])
       AC_CHECK_LIB([Frame],[FrLibIni],,[AC_MSG_ERROR([could not find the frame library])])
       AC_CHECK_HEADERS([FrameL.h],,[AC_MSG_ERROR([could not find the FrameL.h header])])
       AC_DEFINE([USE_FRAMEL],[],[Define to use FrameL])
       FRAME_LIBRARY="FrameL"])

# check for lal libraries and headers
LALSUITE_CHECK_LIB([LAL],[6.12.0],[LALVersion],[lal/LALStdio.h])
LALSUITE_CHECK_LIB([LALSupport],[6.12.0],[XLALCacheImport],[lal/LALCache.h])

# check for inotify
if test "${inotify}" = "true"; then
  AC_CHECK_HEADER([sys/inotify.h],
                  [AC_DEFINE([HAVE_INOTIFY],[1],[Define if you have inotify])],
                  [inotify=false])
fi
INOTIFY_ENABLE_VAL="`eval test x$inotify = xtrue && echo "ENABLED" || echo "DISABLED"`"

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
LALFrame has now been successfully configured:

* using the $FRAME_LIBRARY frame library
* inotify support is $INOTIFY_ENABLE_VAL
* SWIG wrapping module for Octave is $SWIG_OCTAVE_ENABLE_VAL
* SWIG wrapping module for Python is $SWIG_PYTHON_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALFrame,
and run 'make install' to install LALFrame.
==================================================
])

AS_IF([test "${FRAME_LIBRARY}" = "FrameC"],[
       echo "
WARNING: You have selected FrameC as the frame library, support for
FrameC is _not_ complete and is under active development.
"])
