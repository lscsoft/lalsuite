# lal.m4 - lal specific macros
#
# serial 10

AC_DEFUN([LAL_WITH_EXTRA_CPPFLAGS],
[AC_ARG_WITH(
  [extra_cppflags],
  AC_HELP_STRING([--with-extra-cppflags=CPPFLAGS],[additional C preprocessor flags]),
  AS_IF([test -n "${with_extra_cppflags}"],[CPPFLAGS="$CPPFLAGS ${with_extra_cppflags}"]),)
])

AC_DEFUN([LAL_WITH_CFLAGS],
[AC_ARG_WITH(
  [cflags],
  AC_HELP_STRING([--with-cflags=CFLAGS],[C compiler flags]),
  [ if test -n "${with_cflags}"
    then
      CFLAGS="${with_cflags}";
    fi
  ],)
])

AC_DEFUN([LAL_WITH_EXTRA_CFLAGS],
[AC_ARG_WITH(
  [extra_cflags],
  AC_HELP_STRING([--with-extra-cflags=CFLAGS],[additional C compiler flags]),
  [ if test -n "${with_extra_cflags}"
    then
      CFLAGS="$CFLAGS ${with_extra_cflags}";
    fi
  ],)
])

AC_DEFUN([LAL_WITH_EXTRA_LDFLAGS],
[AC_ARG_WITH(
  [extra_ldflags],
  AC_HELP_STRING([--with-extra-ldflags=LDFLAGS],[additional linker flags]),
  [ if test -n "${with_extra_ldflags}"
    then
      LDFLAGS="$LDFLAGS ${with_extra_ldflags}";
    fi
  ],)
])

AC_DEFUN([LAL_WITH_EXTRA_LIBS],
[AC_ARG_WITH(
  [extra_libs],
  AC_HELP_STRING([--with-extra-libs=LIBS],[additional -l and -L linker flags]),
  [ if test -n "${with_extra_libs}"
    then
      LIBS="$LIBS ${with_extra_libs}";
    fi
  ],)
])

AC_DEFUN([LAL_WITH_CC],
[AC_ARG_WITH(
  [cc],
  AC_HELP_STRING([--with-cc=CC],[use the CC C compiler]),
  [ if test -n "${with_cc}"
    then
      CC="${with_cc}";
    fi
  ],)
])

AC_DEFUN([LAL_ENABLE_INTELFFT],
[AC_ARG_ENABLE(
  [intelfft],
  AC_HELP_STRING([--enable-intelfft],[use Intel FFT libraries insted of FFTW [default=no]]),
  [ case "${enableval}" in
      yes) intelfft=true ;;
      no)  intelfft=false ;;
      condor) intelfft=true ; qthread=true ; AC_DEFINE(LAL_QTHREAD, 1, Use fake qthread library for MKL Condor compatibility) ;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-intelfft) ;;
    esac
  ], [ intelfft=false ] )
])

AC_DEFUN([LAL_ENABLE_MACROS],
[AC_ARG_ENABLE(
  [macros],
  AC_HELP_STRING([--enable-macros],[use LAL macros [default=yes]]),
  [ case "${enableval}" in
      yes) ;;
      no) AC_DEFINE(NOLALMACROS, 1, Use functions rather than macros) ;;
      *) AC_MSG_ERROR(bad value for ${enableval} for --enable-debug) ;;
    esac
  ], )
])

AC_DEFUN([LAL_ENABLE_PTHREAD_LOCK],
[AC_ARG_ENABLE(
  [pthread_lock],
  AC_HELP_STRING([--enable-pthread-lock],[use pthread mutex lock for threadsafety [default=no]]),
  [ case "${enableval}" in
      yes) lal_pthread_lock=true; AC_DEFINE(LAL_PTHREAD_LOCK, 1, Use pthread mutex lock for threadsafety) ;;
      no) lal_pthread_lock=false ;;
      *) AC_MSG_ERROR(bad value for ${enableval} for --enable-pthread-lock) ;;
    esac
  ], [ lal_pthread_lock=false ] )
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
