# AX_PYTHON_DEVEL(VERSION)
#
AC_DEFUN([AX_PYTHON_DEVEL],
[if test -n "$1" ; then
  for version in $1 ; do
    AC_MSG_CHECKING([Checking for python$version-config])
    AC_PATH_TOOL([PYTHON_CONFIG], [python$version-config])
    if test -n "$PYTHON_CONFIG" ; then
      AC_MSG_RESULT([yes])
      PYTHON_VERSION=$version
      break
    fi
    AC_MSG_RESULT([no])
  done
else
  AC_MSG_CHECKING([Checking for python-config])
  AC_PATH_TOOL([PYTHON_CONFIG], [python-config])
  if test -n "$PYTHON_CONFIG" ; then
    AC_MSG_RESULT([yes])
    PYTHON_VERSION=$(python -c "import sys;print '.'.join(map(str, sys.version_info@<:@:2@:>@))")
  else
    AC_MSG_RESULT([no])
  fi
fi

if test -z "$PYTHON_CONFIG" ; then
  AC_MSG_ERROR([no python-config found.])
fi

PYTHON_CPPFLAGS="$("$PYTHON_CONFIG" --includes)"
PYTHON_LDFLAGS="$("$PYTHON_CONFIG" --ldflags)"

ac_save_CPPFLAGS="$CPPFLAGS"
ac_save_LIBS="$LIBS"
CPPFLAGS="$LIBS $PYTHON_CPPFLAGS"
LIBS="$LIBS $PYTHON_LIBS"
AC_MSG_CHECKING([Checking if python-config results are accurate])
AC_LANG_PUSH([C])
AC_LINK_IFELSE([
  AC_LANG_PROGRAM([[#include <Python.h>]],
                  [[Py_Initialize();]])
  ],
  [AC_MSG_RESULT([yes])]
  [AC_MSG_RESULT([no])
   AC_MSG_ERROR([$PYTHON_CONFIG output is not usable])])
AC_LANG_POP([C])

CPPFLAGS="$ac_save_CPPFLAGS"
LIBS="$ac_save_LIBS"

AC_SUBST([PYTHON_VERSION])
AC_SUBST([PYTHON_CPPFLAGS])
AC_SUBST([PYTHON_LDFLAGS])

])
