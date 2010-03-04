AC_DEFUN([LAL_ENABLE_GCC_FLAGS],
[AC_ARG_ENABLE([gcc_flags],
  AC_HELP_STRING([--enable-gcc-flags],[turn on strict gcc warning flags (default=yes)]),
  [case "${enableval}" in
     yes) DO_ENABLE_LAL_GCC_FLAGS;;
     no) ;;
     *) DO_ENABLE_LAL_GCC_FLAGS;;
   esac ],
   [ DO_ENABLE_LAL_GCC_FLAGS ] )
])

AC_DEFUN([DO_ENABLE_LAL_GCC_FLAGS],
[
  lal_gcc_flags="-g3 -O4 -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fno-common -Wnested-externs -Wno-format-zero-length"
  case $host_cpu-$host_os in
    *i386-darwin*) lal_gcc_flags="${lal_gcc_flags} -pedantic" ;;
    *x86_64-darwin*) lal_gcc_flags="${lal_gcc_flags} -pedantic" ;;
    *) lal_gcc_flags="${lal_gcc_flags} -pedantic-errors" ;;
  esac
])

AC_DEFUN([LAL_ENABLE_BOINC],
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

AC_DEFUN([LAL_CHECK_BOINC],
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

