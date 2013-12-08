AC_PREREQ([2.63])
AC_INIT([lalmetaio],[1.1.1.1],[lal-discuss@ligo.org])
AC_CONFIG_HEADERS([src/config.h src/LALMetaIOConfig.h])
AC_CONFIG_SRCDIR([src/LIGOLwXML.c])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([ \
  lalmetaio.spec \
  lalmetaio.pc \
  lalmetaio-uninstalled.pc \
  Makefile \
  gnuscripts/Makefile \
  debian/Makefile \
  include/Makefile \
  include/lal/Makefile \
  python/Makefile \
  python/__init__.py \
  src/Makefile \
  swig/Makefile \
  swig/octave/Makefile \
  swig/python/Makefile
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
AC_SUBST([LIBCURRENT],[3])
AC_SUBST([LIBREVISION],[1])
AC_SUBST([LIBAGE],[0])
AC_SUBST([LIBVERSION],[${LIBCURRENT}:${LIBREVISION}:${LIBAGE}])

# nightly build
LALSUITE_ENABLE_NIGHTLY

# debug
LALSUITE_ENABLE_DEBUG

# enable strict gcc flags
LALSUITE_ENABLE_GCC_FLAGS

# enable SWIG wrapping modules
LALSUITE_ENABLE_SWIG

# put version/configure info in config header
lalmetaio_version_major=`echo "$VERSION" | cut -d. -f1`
lalmetaio_version_minor=`echo "$VERSION" | cut -d. -f2`
lalmetaio_version_micro=`echo "$VERSION" | cut -d. -f3`
lalmetaio_version_devel=`echo "$VERSION" | cut -d. -f4-`
test -z "$lalmetaio_version_micro" && lalmetaio_version_micro=0
test -z "$lalmetaio_version_devel" && lalmetaio_version_devel=0
lalmetaio_configure_date=`date +"%Y-%m-%dT%H:%M:%S%z"`
AC_DEFINE_UNQUOTED([LALMETAIO_VERSION],["$VERSION"],[LALMetaIO Version])
AC_DEFINE_UNQUOTED([LALMETAIO_VERSION_MAJOR],[$lalmetaio_version_major],
                   [LALMetaIO Version Major Number])
AC_DEFINE_UNQUOTED([LALMETAIO_VERSION_MINOR],[$lalmetaio_version_minor],
                   [LALMetaIO Version Minor Number])
AC_DEFINE_UNQUOTED([LALMETAIO_VERSION_MICRO],[$lalmetaio_version_micro],
                   [LALMetaIO Version Micro Number])
AC_DEFINE_UNQUOTED([LALMETAIO_VERSION_DEVEL],[$lalmetaio_version_devel],
                   [LALMetaIO Version Devel Number])
AC_DEFINE_UNQUOTED([LALMETAIO_CONFIGURE_ARGS],["$ac_configure_args"],
                   [LALMetaIO Configure Arguments])
AC_DEFINE_UNQUOTED([LALMETAIO_CONFIGURE_DATE],["$lalmetaio_configure_date"],
                   [LALMetaIO Configure Date])

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
AM_PATH_PYTHON([2.6],,[AC_MSG_ERROR([Python-2.6, or higher, is required])])

# check for supported mac os x version
LALSUITE_OSX_VERSION_CHECK

# link tests using libtool
LALSUITE_USE_LIBTOOL

# check for header files
AC_HEADER_STDC
AC_CHECK_HEADERS([zlib.h],,[AC_MSG_ERROR([could not find zlib.h])])

# check for required libraries
AC_CHECK_LIB([m],[main],,[AC_MSG_ERROR([could not find the math library])])
AC_CHECK_LIB([z],[compress],,[AC_MSG_ERROR([could not find the zlib library])])

# metaio
PKG_CHECK_MODULES([METAIO],[libmetaio],[true],[false])
LALSUITE_ADD_FLAGS([C],[${METAIO_CFLAGS}],[${METAIO_LIBS}])
AC_CHECK_LIB([metaio],[MetaioOpen],,[AC_MSG_ERROR([could not find the metaio library])])

# check for libmetaio headers
AC_CHECK_HEADERS([metaio.h],,[AC_MSG_ERROR([could not find the metaio.h header])])

# check for lal libraries and headers
LALSUITE_CHECK_LIB([LAL],[6.11.0],[LALVersion],[lal/LALStdio.h])
LALSUITE_CHECK_LIB([LALSupport],[6.11.0],[LALOpenDataFile],[lal/FileIO.h])

# configure SWIG wrapping modules
LALSUITE_USE_SWIG([XLAL LAL])

# add gcc specific flags
LALSUITE_ADD_GCC_FLAGS

# force all python modules into the exec dir to prevent mixed C and
# pure-python library code from being split across lib and lib64 on
# RedHat-derived systems.
AC_SUBST([pkgpythondir], ["\${pkgpyexecdir}"])

LAL_HEADER_LINK=${ac_aux_dir}/header-link.mk
AC_SUBST_FILE(LAL_HEADER_LINK)

# output
AC_OUTPUT
echo "
==================================================
LALMetaIO has now been successfully configured:

* SWIG wrapping module for Octave is $SWIG_OCTAVE_ENABLE_VAL
* SWIG wrapping module for Python is $SWIG_PYTHON_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALMetaIO,
and run 'make install' to install LALMetaIO.
==================================================
"
