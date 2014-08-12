AC_PREREQ([2.63])
AC_INIT([LALFrame],[1.1.1.1],[lal-discuss@ligo.org])
AC_CONFIG_HEADERS([src/config.h src/LALFrameConfig.h])
AC_CONFIG_SRCDIR([src/LALFrameIO.c])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([ \
  lalframe.spec \
  lalframe.pc \
  lalframe-uninstalled.pc \
  Makefile \
  debian/Makefile \
  include/Makefile \
  include/lal/Makefile \
  src/LALFrameVCSInfo.c \
  src/LALFrameVCSInfo.h \
  src/Makefile \
  python/Makefile \
  python/lalframe/Makefile \
  python/lalframe/__init__.py \
  python/lalframe/utils/Makefile \
  swig/Makefile \
  test/Makefile
])
AM_INIT_AUTOMAKE([1.11 foreign subdir-objects color-tests parallel-tests])

AC_CANONICAL_HOST

LALSUITE_ARG_VAR

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

# enable check for FrameC
LALFRAME_ENABLE_FRAMEC

# enable check for inotify
LALFRAME_ENABLE_INOTIFY

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

# check for pythread, needed for low latency data test codes
AX_PTHREAD([lalframe_pthread=true],[lalframe_pthread=false])
AM_CONDITIONAL([PTHREAD],[test x$lalframe_pthread = xtrue])

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

# checks for library functions
AC_CHECK_FUNCS([gmtime_r localtime_r])

# check for framec or libframe libraries and headers
FRAMEC_AVAILABLE="no"
FRAMEL_AVAILABLE="yes"
PKG_PROG_PKG_CONFIG
AS_IF([test "${framec}" = "true"],[
	FRAMEC_AVAILABLE="yes"
	PKG_CHECK_MODULES([FRAMECPPC],[framecppc >= 1.19.32],[true],[false])
	LALSUITE_ADD_FLAGS([C],[${FRAMECPPC_CFLAGS}],[${FRAMECPPC_LIBS}])
	AC_CHECK_LIB([framecppc],[FrameCFileOpen],,[FRAMEC_AVAILABLE="no"])
	AC_CHECK_HEADERS([framecppc/FrameC.h],,[FRAMEC_AVAILABLE="no"])])
PKG_CHECK_MODULES([LIBFRAME],[libframe],[true],[false])
LALSUITE_ADD_FLAGS([C],[${LIBFRAME_CFLAGS}],[${LIBFRAME_LIBS}])
AC_CHECK_LIB([Frame],[FrLibIni],,[FRAMEL_AVAILABLE="no"])
AC_CHECK_HEADERS([FrameL.h],,[FRAMEL_AVAILABLE="no"])

# make sure a frame library was found
AS_IF([test "${FRAMEC_AVAILABLE}" = "no" -a "${FRAMEL_AVAILABLE}" = "no"], [
	AC_MSG_ERROR([no frame library available])
])
AM_CONDITIONAL([FRAMEC],[test "${FRAMEC_AVAILABLE}" = "yes"])
AM_CONDITIONAL([FRAMEL],[test "${FRAMEL_AVAILABLE}" = "yes"])

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
LALSUITE_USE_SWIG

# add gcc specific flags
LALSUITE_ADD_GCC_FLAGS

# link any executables against LALFrame library
AC_SUBST([LDADD], ["\$(top_builddir)/src/liblalframe.la"])

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
LALFrame has now been successfully configured:

* Python support is $PYTHON_ENABLE_VAL
* using the $FRAME_LIBRARY frame library
* FrameC availability... ${FRAMEC_AVAILABLE}
* FrameL availability... ${FRAMEL_AVAILABLE}
* inotify support is $INOTIFY_ENABLE_VAL
* SWIG binding generation is $SWIG_GENERATE_ENABLE_VAL
* SWIG bindings for Octave are $SWIG_BUILD_OCTAVE_ENABLE_VAL
* SWIG bindings for Python are $SWIG_BUILD_PYTHON_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALFrame,
and run 'make install' to install LALFrame.
==================================================
])
