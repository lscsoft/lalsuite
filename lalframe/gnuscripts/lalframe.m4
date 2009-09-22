dnl lalxml.m4

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
    *) lal_gcc_flags="${lal_gcc_flags} -pedantic-errors" ;;
  esac
])
