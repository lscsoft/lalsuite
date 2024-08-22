# lal.m4 - lal specific macros
#
# serial 25

AC_DEFUN([LAL_WITH_DEFAULT_DEBUG_LEVEL],[
  AC_ARG_WITH(
    [default_debug_level],
    AS_HELP_STRING([--with-default-debug-level],[set default value of lalDebugLevel [default=1, i.e. print error messages]]),
    [AS_IF([test "x`expr "X${withval}" : ["^\(X[0-9][0-9]*$\)"]`" != "xX${withval}"],[
      AC_MSG_ERROR([bad integer value '${withval}' for --with-default-debug-level])
    ])],
    [with_default_debug_level=1]
  )
  AC_DEFINE_UNQUOTED([LAL_DEFAULT_DEBUG_LEVEL],[${with_default_debug_level}],[Set default value of lalDebugLevel])
])

AC_DEFUN([LAL_ENABLE_FFTW3_MEMALIGN],
[AC_ARG_ENABLE(
  [fftw3_memalign],
  AS_HELP_STRING([--enable-fftw3-memalign],[use aligned memory optimizations with fftw3 [default=no]]),
  AS_CASE(["${enableval}"],
    [yes],[fftw3_memalign=true],
    [no],[fftw3_memalign=false],
    AC_MSG_ERROR([bad value for ${enableval} for --enable-fftw3-memalign])
  ),[fftw3_memalign=false])
])

AC_DEFUN([LAL_ENABLE_INTELFFT],
[AC_ARG_ENABLE(
  [intelfft],
  AS_HELP_STRING([--enable-intelfft],[use Intel FFT libraries insted of FFTW [default=no]]),
  AS_CASE(["${enableval}"],
    [yes],[intelfft=true],
    [no],[intelfft=false],
    AC_MSG_ERROR([bad value for ${enableval} for --enable-intelfft])
  ),[intelfft=false])
])

AC_DEFUN([LAL_ENABLE_MEMORY_FUNCTIONS],
[AC_ARG_ENABLE(
  [memory-functions],
  AS_HELP_STRING([--enable-memory-functions],[use LAL memory functions [default=yes]]),
  AS_CASE(["${enableval}"],
    [yes],,
    [no],[AC_DEFINE([LAL_MEMORY_FUNCTIONS_DISABLED],[1],[Disable LAL memory functions])],
    AC_MSG_ERROR([bad value for ${enableval} for --enable-memory-functions])
  ),)
])

AC_DEFUN([LAL_ENABLE_PTHREAD_LOCK], [
  AC_ARG_ENABLE([pthread_lock],
    AS_HELP_STRING([--enable-pthread-lock],[use pthread mutex lock for threadsafety @<:@default=yes@:>@]),
    AS_CASE(["${enableval}"],
      [yes],[lal_pthread_lock=true],
      [no],[lal_pthread_lock=false],
      AC_MSG_ERROR([bad value for ${enableval} for --enable-pthread-lock])
    ),
    [lal_pthread_lock=true]
  )
  AS_IF([test "x$lal_pthread_lock" = "xtrue"], [
    # check for pthreads
    AX_PTHREAD([
      LALSUITE_ADD_FLAGS([C],[${PTHREAD_CFLAGS}],[${PTHREAD_LIBS}])
    ],[
      AC_MSG_ERROR([do not know how to compile posix threads])
    ])
    AC_DEFINE([LAL_PTHREAD_LOCK],[1],[Use pthread mutex lock for threadsafety])
  ])
  AC_SUBST([PTHREAD_CFLAGS])
  AC_SUBST([PTHREAD_LIBS])
])

AC_DEFUN([LAL_INTEL_FFT_LIBS_MSG_ERROR],
[echo "**************************************************************"
 echo "* The --enable-static and --enable-shared options are        *"
 echo "* mutually exclusive with Intel FFT libraries.               *"
 echo "*                                                            *"
 echo "* Please reconfigure with:                                   *"
 echo "*                                                            *"
 echo "*   --enable-static --disable-shared                         *"
 echo "*                                                            *"
 echo "* for static libraries or                                    *"
 echo "*                                                            *"
 echo "*   --disable-static --enable-shared                         *"
 echo "*                                                            *"
 echo "* for shared libraries.                                      *"
 echo "*                                                            *"
 echo "* Please see the instructions in the file INSTALL.           *"
 echo "**************************************************************"
AC_MSG_ERROR([Intel FFT must use either static or shared libraries])
])
