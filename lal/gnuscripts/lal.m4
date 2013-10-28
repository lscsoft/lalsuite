# lal.m4 - lal specific macros
#
# serial 17

AC_DEFUN([LAL_ENABLE_INTELFFT],
[AC_ARG_ENABLE(
  [intelfft],
  AC_HELP_STRING([--enable-intelfft],[use Intel FFT libraries insted of FFTW [default=no]]),
  AS_CASE(["${enableval}"],
    [yes],[intelfft=true],
    [no],[intelfft=false],
    [condor],[intelfft=true; qthread=tru; AC_DEFINE([LALQTHREAD],[1],[Use fake qthread library for MKL Condor compatibility])],
    AC_MSG_ERROR([bad value for ${enableval} for --enable-intelfft])
  ),[intelfft=false])
])

AC_DEFUN([LAL_ENABLE_MACROS],
[AC_ARG_ENABLE(
  [macros],
  AC_HELP_STRING([--enable-macros],[use LAL macros [default=yes]]),
  AS_CASE(["${enableval}"],
    [yes],,
    [no],AC_DEFINE([NOLALMACROS],[1],[Use functions rather than macros]),
    AC_MSG_ERROR([bad value for ${enableval} for --enable-macros])
  ),)
])

AC_DEFUN([LAL_ENABLE_PTHREAD_LOCK],
[AC_ARG_ENABLE(
  [pthread_lock],
  AC_HELP_STRING([--enable-pthread-lock],[use pthread mutex lock for threadsafety [default=no]]),
  AS_CASE(["${enableval}"],
    [yes],[lal_pthread_lock=true; AC_DEFINE([LAL_PTHREAD_LOCK],[1],[Use pthread mutex lock for threadsafety])],
    [no],[lal_pthread_lock=false],
    AC_MSG_ERROR([bad value for ${enableval} for --enable-pthread-lock])
  ),[lal_pthread_lock=false])
])

AC_DEFUN([LAL_INTEL_MKL_QTHREAD_WARNING],
[echo "**************************************************************"
 echo "* LAL will be linked against the fake POSIX thread library!  *"
 echo "*                                                            *"
 echo "* This build of LAL will not be thread safe and cannot be    *"
 echo "* linked against the system pthread library.                 *"
 echo "*                                                            *"
 echo "* The environment variables                                  *"
 echo "*                                                            *"
 echo "*    MKL_SERIAL=YES                                          *"
 echo "*    KMP_LIBRARY=serial                                      *"
 echo "*                                                            *"
 echo "* must be set before running executables linked against this *"
 echo "* build of LAL.                                              *"
 echo "*                                                            *"
 echo "* Please see the documention of the FFT package for details. *"
 echo "**************************************************************"
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

AC_DEFUN([LAL_CHECK_GSL_VERSION],
[
  lal_min_gsl_version=ifelse([$1], ,1.0,$1)
  AC_MSG_CHECKING(for GSL version >= $lal_min_gsl_version)
  AC_TRY_RUN([
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_version.h>
int main(void)
{
  int required_major, required_minor;
  int major, minor;
  char required_version[] = "$lal_min_gsl_version";
  char version[] = GSL_VERSION;
  if ( strcmp(GSL_VERSION, gsl_version) ) {
    printf("error\n*** mismatch between header and library versions of GSL\n" );
    printf("\n*** header  has version %s\n", GSL_VERSION);
    printf("\n*** library has version %s\n", gsl_version);
    exit(1);
  }
  sscanf(required_version, "%d.%d", &required_major, &required_minor);
  sscanf(version, "%d.%d", &major, &minor);
  if ( major < required_major || (major == required_major && minor < required_minor) ) {
    printf("no\n*** found version %s of GSL but minimum version is %d.%d\n", GSL_VERSION, required_major, required_minor );
    exit(1);
  }
  return 0;
}
  ], [AC_MSG_RESULT(yes)], [AC_MSG_ERROR(could not find required version of GSL)], [echo $ac_n "cross compiling; assumed OK... $ac_c"])
])
