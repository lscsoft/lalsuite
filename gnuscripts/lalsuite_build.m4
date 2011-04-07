# lalsuite_build.m4 - top level build macros
#
# serial 15

AC_DEFUN([LALSUITE_USE_LIBTOOL],
[## $0: Generate a libtool script for use in configure tests
AC_PROVIDE_IFELSE([LT_INIT], ,
                  [m4_fatal([$0: requires libtool])])[]dnl
LT_OUTPUT
m4_append([AC_LANG(C)],
[ac_link="./libtool --mode=link --tag=CC $ac_link"
])[]dnl
AC_PROVIDE_IFELSE([AC_PROG_CXX],
[m4_append([AC_LANG(C++)],
[ac_link="./libtool --mode=link --tag=CXX $ac_link"
])])[]dnl
AC_LANG(_AC_LANG)[]dnl
]) # LALSUITE_USE_LIBTOOL

AC_DEFUN([LALSUITE_ARG_VAR],[
  AC_ARG_VAR(LALSUITE_BUILD,[Set if part of lalsuite build])
  AC_ARG_VAR(LALSUITE_TOP_SRCDIR,[Set to top source directory of lalsuite])
])

AC_DEFUN([LALSUITE_ENABLE_MODULE],[
AM_CONDITIONAL([$1],[test x$$2 = xtrue])
eval $1_ENABLE_VAL="`eval test "$$2" = "true" && echo "ENABLED" || echo "DISABLED"`"
])

AC_DEFUN([LALSUITE_CHECK_LIB],[
m4_pushdef([lowercase],translit([[$1]], [A-Z], [a-z]))
m4_pushdef([uppercase],translit([[$1]], [a-z], [A-Z]))
PKG_CHECK_MODULES(uppercase,[lowercase >= $2],[lowercase="true"],[lowercase="false"])
if test "$lowercase" = "true"; then
  CPPFLAGS="$CPPFLAGS $[]uppercase[]_CFLAGS"
  LIBS="$LIBS $[]uppercase[]_LIBS"
  if test "$LALSUITE_BUILD" = "true"; then
    AC_DEFINE([HAVE_LIB]uppercase,[1],[Define to 1 if you have the $1 library])
    lowercase="true"
  else
    AC_CHECK_LIB(lowercase,[$3],[lowercase="true"],[AC_MSG_ERROR([could not find the $1 library])])
    AC_CHECK_HEADERS([$4],,[AC_MSG_ERROR([could not find the $4 header])])
    if test "$1" != "LALSupport"; then
      LALSUITE_HEADER_LIBRARY_MISMATCH_CHECK([$1])
    fi
    AC_DEFINE([HAVE_LIB]uppercase,[1],[Define to 1 if you have the $1 library])
  fi
else
  AC_MSG_ERROR([could not find the $1 library])
fi
LALSUITE_ENABLE_MODULE(uppercase,lowercase)
m4_popdef([lowercase])
m4_popdef([uppercase])
])

AC_DEFUN([LALSUITE_CHECK_OPT_LIB],[
m4_pushdef([lowercase],translit([[$1]], [A-Z], [a-z]))
m4_pushdef([uppercase],translit([[$1]], [a-z], [A-Z]))
if test "$lowercase" = "true"; then
  PKG_CHECK_MODULES(uppercase,[lowercase >= $2],[lowercase="true"],[lowercase="false"])
  if test "$lowercase" = "true"; then
    if test "$LALSUITE_BUILD" = "true"; then
      AC_DEFINE([HAVE_LIB]uppercase,[1],[Define to 1 if you have the $1 library])
      lowercase="true"
      CPPFLAGS="$CPPFLAGS $[]uppercase[]_CFLAGS"
      LIBS="$LIBS $[]uppercase[]_LIBS"
    else
      CPPFLAGS="$CPPFLAGS $[]uppercase[]_CFLAGS"
      LIBS="$LIBS $[]uppercase[]_LIBS"
      AC_CHECK_LIB(lowercase,[$3],[lowercase="true"],[lowercase=false
        AC_MSG_WARN([could not find the $1 library])])
      if test "$lowercase" = true; then
        AC_CHECK_HEADERS([$4],,[lowercase=false])
        if test "$lowercase" = true; then
          if test "$1" != "LALSupport"; then
            LALSUITE_HEADER_LIBRARY_MISMATCH_CHECK([$1])
          fi
          if test "$lowercase" = true; then
            AC_DEFINE([HAVE_LIB]uppercase,[1],[Define to 1 if you have the $1 library])
          fi
        fi
      fi
    fi
  fi
fi
LALSUITE_ENABLE_MODULE(uppercase,lowercase)
m4_popdef([lowercase])
m4_popdef([uppercase])
])

AC_DEFUN([LALSUITE_HEADER_LIBRARY_MISMATCH_CHECK],[
AC_MSG_CHECKING([whether $1 headers match the library])
lib_structure=`echo $1 | sed 's/LAL/lal/'`VCSInfo
header_structure=`echo $1 | sed 's/LAL/lal/'`HeaderVCSInfo
AC_RUN_IFELSE(
  [AC_LANG_SOURCE([[
#include <string.h>
#include <stdlib.h>
#include <lal/$1VCSInfo.h>
int main(void) { exit(XLALVCSInfoCompare(&$lib_structure, &$header_structure) ? 1 : 0); }
  ]])],
  [
    AC_MSG_RESULT(yes)
  ],
  [
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([Your $1 headers do not match your
library. Check config.log for details.
])
  ],
  [
    AC_MSG_WARN([cross compiling: not checking])
  ]
)
])

AC_DEFUN([LALSUITE_ENABLE_NIGHTLY],
[AC_ARG_ENABLE(
  [nightly],
  AC_HELP_STRING([--enable-nightly],[nightly build [default=no]]),
  [ case "${enableval}" in
      yes) NIGHTLY_VERSION=`date +"%Y%m%d"`
           VERSION="${VERSION}.${NIGHTLY_VERSION}" ;;
      no) NIGHTLY_VERSION="";;
      *) NIGHTLY_VERSION="${enableval}"
         VERSION="${VERSION}.${NIGHTLY_VERSION}" ;;
      esac ],
  [ NIGHTLY_VERSION="" ] )
  AC_SUBST(NIGHTLY_VERSION)
])

AC_DEFUN([LALSUITE_ENABLE_LALFRAME],
[AC_ARG_ENABLE(
  [lalframe],
  AC_HELP_STRING([--enable-lalframe],[compile code that requires lalframe library [default=yes]]),
  [ case "${enableval}" in
      yes) lalframe=true;;
      no) lalframe=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalframe) ;;
    esac
  ], [ lalframe=true ] )
if test "$frame" = "false"; then
  lalframe=false
fi
])

AC_DEFUN([LALSUITE_ENABLE_LALMETAIO],
[AC_ARG_ENABLE(
  [lalmetaio],
  AC_HELP_STRING([--enable-lalmetaio],[compile code that requires lalmetaio library [default=yes]]),
  [ case "${enableval}" in
      yes) lalmetaio=true;;
      no) lalmetaio=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalmetaio) ;;
    esac
  ], [ lalmetaio=true ] )
if test "$metaio" = "false"; then
  lalmetaio=false
fi
])

AC_DEFUN([LALSUITE_ENABLE_LALXML],
[AC_ARG_ENABLE(
  [lalxml],
  AC_HELP_STRING([--enable-lalxml],[compile code that requires lalxml library [default=no]]),
  [ case "${enableval}" in
      yes) lalxml=true;;
      no) lalxml=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalxml) ;;
    esac
  ], [ lalxml=false ] )
])

AC_DEFUN([LALSUITE_ENABLE_LALBURST],
[AC_ARG_ENABLE(
  [lalburst],
  AC_HELP_STRING([--enable-lalburst],[compile code that requires lalburst library [default=yes]]),
  [ case "${enableval}" in
      yes) lalburst=true;;
      no) lalburst=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalburst) ;;
    esac
  ], [ lalburst=true ] )
if test "$lalmetaio" = "false"; then
  lalburst=false
fi])

AC_DEFUN([LALSUITE_ENABLE_LALINSPIRAL],
[AC_ARG_ENABLE(
  [lalinspiral],
  AC_HELP_STRING([--enable-lalinspiral],[compile code that requires lalinspiral library [default=yes]]),
  [ case "${enableval}" in
      yes) lalinspiral=true;;
      no) lalinspiral=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalinspiral) ;;
    esac
  ], [ lalinspiral=true ] )
if test "$lalmetaio" = "false"; then
  lalinspiral=false
fi
])

AC_DEFUN([LALSUITE_ENABLE_LALPULSAR],
[AC_ARG_ENABLE(
  [lalpulsar],
  AC_HELP_STRING([--enable-lalpulsar],[compile code that requires lalpulsar library [default=yes]]),
  [ case "${enableval}" in
      yes) lalpulsar=true;;
      no) lalpulsar=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalpulsar) ;;
    esac
  ], [ lalpulsar=true ] )
])

AC_DEFUN([LALSUITE_ENABLE_LALSTOCHASTIC],
[AC_ARG_ENABLE(
  [lalstochastic],
  AC_HELP_STRING([--enable-lalstochastic],[compile code that requires lalstochastic library [default=yes]]),
  [ case "${enableval}" in
      yes) lalstochastic=true;;
      no) lalstochastic=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalstochastic) ;;
    esac
  ], [ lalstochastic=true ] )
if test "$lalmetaio" = "false"; then
  lalstochastic=false
fi
])

AC_DEFUN([LALSUITE_ENABLE_BOINC],
[AC_ARG_ENABLE(
  [boinc],
  AC_HELP_STRING([--enable-boinc],[enable BOINC support [default=no]]),
  [ case "${enableval}" in
      yes) boinc=true;;
      no) boinc=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-boinc);;
    esac
  ], [ boinc=false ] )
AC_ARG_VAR([BOINC_PREFIX],[BOINC installation directory (optional)])
])

AC_DEFUN([LALSUITE_CHECK_BOINC],
[AC_MSG_CHECKING([whether LAL has been compiled with BOINC support])
AC_TRY_RUN([
#include <lal/LALConfig.h>
#ifdef LAL_BOINC_ENABLED
int main( void ) { return 0; }
#else
int main( void ) { return 1; }
#endif
],
AC_MSG_RESULT([yes])
[boinc=true],
AC_MSG_RESULT([no])
[boinc=false],
AC_MSG_RESULT([unknown])
[boinc=false])
])

AC_DEFUN([LALSUITE_WITH_CUDA],
[AC_ARG_WITH(
  [cuda],
  AC_HELP_STRING([--with-cuda=PATH],[specify location of CUDA [/opt/cuda]]),
  [ case "$with_cuda" in
    no)
      cuda=false
      ;;
    yes)
      if test "x$build_os" != "xlinux"; then
        AC_MSG_ERROR([CUDA not supported on this platform])
      else
        AC_MSG_WARN([No path for CUDA specifed, using /opt/cuda])
        cuda=true
        if test "x$build_cpu" = "xx86_64"; then
          CLIBS="lib64"
        else
          CLIBS="lib"
        fi
        CUDA_LIBS="-L/opt/cuda/$CLIBS -lcufft -lcudart"
        CUDA_CFLAGS="-I/opt/cuda/include"
        LIBS="$LIBS $CUDA_LIBS"
        CFLAGS="$CFLAGS $CUDA_CFLAGS"
        AC_SUBST(CUDA_LIBS)
        AC_SUBST(CUDA_CFLAGS)
      fi
      ;;
    *)
      AC_MSG_NOTICE([Using ${with_cuda} as CUDA path])
      cuda=true
      if test "x$build_cpu" = "xx86_64"; then
        CLIBS="lib64"
      else
        CLIBS="lib"
      fi
      CUDA_LIBS="-L${with_cuda}/$CLIBS -lcufft -lcudart"
      CUDA_CFLAGS="-I${with_cuda}/include"
      LIBS="$LIBS $CUDA_LIBS"
      CFLAGS="$CFLAGS $CUDA_CFLAGS"
      AC_SUBST(CUDA_LIBS)
      AC_SUBST(CUDA_CFLAGS)
      ;;
    esac
  ], [ cuda=false ])
  LALSUITE_ENABLE_MODULE([CUDA],[cuda])
])

AC_DEFUN([LALSUITE_ENABLE_OSX_VERSION_CHECK],
[AC_ARG_ENABLE(
  [osx_version_check],
  AC_HELP_STRING([--enable-osx-version-check][disable OS X version check [default=yes]]),
  [ case "${enableval}" in
      yes) osx_version_check=true;;
      no) osx_version_check=false;;
      *) AC_MSG_ERROR([bad value ${enableval} for --enable-osx-version-check]);;
    esac
  ], [ osx_version_check=true ] )
])

AC_DEFUN([LALSUITE_OSX_VERSION_CHECK],
[
LALSUITE_ENABLE_OSX_VERSION_CHECK
if test "x${osx_version_check}" = "xtrue"; then
  if test "x$build_vendor" = "xapple"; then
    AC_CHECK_PROGS([SW_VERS],[sw_vers])
    if test "x$SW_VERS" != "x"; then
      AC_MSG_CHECKING([Mac OS X version])
      MACOSX_VERSION=`$SW_VERS -productVersion`
      AC_MSG_RESULT([$MACOSX_VERSION])
    fi
    case "$MACOSX_VERSION" in
      10.0*|10.1*|10.2*|10.3*)
        AC_MSG_ERROR([This version of Mac OS X is not supported])
        ;;
      10.4*|10.5*|10.6*)
        # supported version
        ;;
      *)
        AC_MSG_WARN([Unknown Mac OS X version])
        ;;
    esac
  fi
fi
])
