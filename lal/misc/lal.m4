dnl lal.m4
dnl $Id$

AC_DEFUN([LAL_WITH_GCC_FLAGS],
[AC_ARG_WITH(
        [gcc_flags],   
        [  --with-gcc-flags        turn on strict gcc warning flags],
        [ if test -n "${with_gcc_flags}"
          then
            lal_gcc_flags="-g3 -O4 -ansi -pedantic -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -fno-common -Wnested-externs -D__NO_STRING_INLINES -Wno-long-long"
          else
            lal_gcc_flags=""
          fi
        ], [ lal_gcc_flags="" ] )
])

AC_DEFUN([LAL_WITH_EXTRA_CPPFLAGS],
[AC_ARG_WITH(
	[extra_cppflags], 
        [  --with-extra-cppflags=CPPFLAGS  additional C preprocessor flags],
	[ if test -n "${with_extra_cppflags}"
	  then
	    CPPFLAGS="$CPPFLAGS ${with_extra_cppflags}";
	  fi
	],)
])

AC_DEFUN([LAL_WITH_CFLAGS],
[AC_ARG_WITH(
	[cflags], 
        [  --with-cflags=CFLAGS        C compiler flags],
	[ if test -n "${with_cflags}"
	  then
	    CFLAGS="${with_cflags}";
	  fi
	],)
])

AC_DEFUN([LAL_WITH_EXTRA_CFLAGS],
[AC_ARG_WITH(
	[extra_cflags], 
        [  --with-extra-cflags=CFLAGS  additional C compiler flags],
	[ if test -n "${with_extra_cflags}"
	  then
	    CFLAGS="$CFLAGS ${with_extra_cflags}";
	  fi
	],)
])

AC_DEFUN([LAL_WITH_EXTRA_LDFLAGS],
[AC_ARG_WITH(
	[extra_ldflags], 
        [  --with-extra-ldflags=LDFLAGS  additional linker flags],
	[ if test -n "${with_extra_ldflags}"
	  then
	    LDFLAGS="$LDFLAGS ${with_extra_ldflags}";
	  fi
	],)
])

AC_DEFUN([LAL_WITH_EXTRA_LIBS],
[AC_ARG_WITH(
	[extra_libs], 
        [  --with-extra-libs=LIBS  additional -l and -L linker flags],
	[ if test -n "${with_extra_libs}"
	  then
	    LIBS="$LIBS ${with_extra_libs}";
	  fi
	],)
])

AC_DEFUN([LAL_WITH_MPICC],
[AC_ARG_WITH(
        [mpicc], 
        [  --with-mpicc=MPICC      use the MPICC C compiler for MPI code],
        [ if test -n "${with_mpicc}"
          then
            MPICC="${with_mpicc}";
          fi
        ],)
])

AC_DEFUN([LAL_WITH_CC],
[AC_ARG_WITH(
        [cc], 
        [  --with-cc=CC            use the CC C compiler],
        [ if test -n "${with_cc}"
          then
            CC="${with_cc}";
          fi
        ],)
])


AC_DEFUN([LAL_ENABLE_FRAME],
[AC_ARG_ENABLE(
        [frame],
        [  --enable-frame          compile code that requires Frame library [default=yes] ],
        [ case "${enableval}" in
            yes) frame=true;;
            no)  frame=false ;;
            *) AC_MSG_ERROR(bad value ${enableval} for --enable-frame) ;;
          esac
        ], [ frame=true ] )
])

AC_DEFUN([LAL_ENABLE_MPI],
[AC_ARG_ENABLE(
        [mpi],
        [  --enable-mpi            compile code that requires MPI [default=no] ],
        [ case "${enableval}" in
            yes) mpi=true;;
            no)  mpi=false ;;
            *) AC_MSG_ERROR(bad value ${enableval} for --enable-mpi) ;;
          esac
        ], [ mpi=false ] )
])

AC_DEFUN([LAL_ENABLE_METAIO],
[AC_ARG_ENABLE(
  [metaio],
  [  --enable-metaio         compile code that requires metaio/dataflow library [default=yes] ],
  [ case "${enableval}" in
      yes) metaio=true;;
      no)  metaio=false ;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-frame) ;;
    esac
  ], [ metaio=true ] )
])

AC_DEFUN([LAL_ENABLE_INTELFFT],
[AC_ARG_ENABLE(
        [intelfft],
        [  --enable-intelfft       use Intel FFT libraries insted of FFTW [default=no] ],
        [ case "${enableval}" in
            yes) intelfft=true;;
            no)  intelfft=false ;;
            condor) intelfft=true ; qthread=true ; AC_DEFINE(LAL_QTHREAD, 1, Use fake qthread library for MKL Condor compatibility) ;;
            *) AC_MSG_ERROR(bad value ${enableval} for --enable-intelfft) ;;
          esac
        ], [ intelfft=false ] )
])

AC_DEFUN([LAL_ENABLE_DEBUG],
[AC_ARG_ENABLE(
        [debug],
        [  --enable-debug          include standard LAL debugging code [default=yes] ],
        [ case "${enableval}" in
            yes) ;;
            no)  AC_DEFINE(LAL_NDEBUG, 1, Suppress debugging code) ;;
            *) AC_MSG_ERROR(bad value for ${enableval} for --enable-debug) ;;
          esac
        ],)
])

AC_DEFUN([LAL_ENABLE_MACROS],
[AC_ARG_ENABLE(
        [macros],
        [  --enable-macros         use LAL macros [default=yes] ],
        [ case "${enableval}" in
            yes) ;;
            no)  AC_DEFINE(NOLALMACROS, 1, Use functions rather than macros) ;;
            *) AC_MSG_ERROR(bad value for ${enableval} for --enable-debug) ;;
          esac
        ],)
])

AC_DEFUN([LAL_ENABLE_NIGHTLY],
[AC_ARG_ENABLE(
        [nightly],
        [  --enable-nightly        nightly build [default=no] ],
        [ case "${enableval}" in
            yes) NIGHTLY_VERSION=`date +"%Y_%m_%d"`
                 VERSION="${VERSION}_${NIGHTLY_VERSION}";;
            no)  NIGHTLY_VERSION="";;
            *)   NIGHTLY_VERSION="${enableval}"
                 VERSION="${VERSION}_${NIGHTLY_VERSION}";;
          esac ],
        [ NIGHTLY_VERSION="" ] )
 AC_SUBST(NIGHTLY_VERSION)
])
                

AC_DEFUN([LAL_ENABLE_PTHREAD_LOCK],
[AC_ARG_ENABLE(
        [pthread_lock],
        [  --enable-pthread-lock   use pthread mutex lock for threadsafety [default=no] ],
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


AC_DEFUN([LAL_CHECK_MPI_FLAGS],
[ AC_MSG_CHECKING([for mpicc flags])
  SHOWARG=""
  MPI_CPPFLAGS=""
  MPI_CFLAGS=""
  MPI_LDFLAGS=""
  MPI_TYPE=NONE
  if (($MPICC -compile_info 1>/dev/null 2>/dev/null) && ($MPICC -link_info 1>/dev/null 2>/dev/null)) ; then
    MPI_TYPE=MPICH
    for mpiarg in `$MPICC -compile_info` ; do
      case $mpiarg in
        -D*) MPI_CPPFLAGS="$MPI_CPPFLAGS $mpiarg" ;;
        -I*) MPI_CPPFLAGS="$MPI_CPPFLAGS $mpiarg" ;;
      esac
    done
    for mpiarg in `$MPICC -link_info` ; do
      case $mpiarg in
        -L*) MPI_LDFLAGS="$MPI_LDFLAGS $mpiarg" ;;
        -l*) MPI_LDFLAGS="$MPI_LDFLAGS $mpiarg" ;;
      esac
    done
  else
    if ($MPICC -show 1>/dev/null 2>/dev/null) ; then
      MPI_TYPE=MPICH
      SHOWARG="-show"
    elif ($MPICC -showme 1>/dev/null 2>/dev/null) ; then
      MPI_TYPE=LAM
      SHOWARG="-showme"
    else
      AC_MSG_WARN([couldn't determine mpi compile flags])
    fi
    if test -n "$SHOWARG" ; then
      for mpiarg in `$MPICC $SHOWARG` ; do
        case $mpiarg in
          -D*) MPI_CPPFLAGS="$MPI_CPPFLAGS $mpiarg" ;;
          -I*) MPI_CPPFLAGS="$MPI_CPPFLAGS $mpiarg" ;;
          -L*) MPI_LDFLAGS="$MPI_LDFLAGS $mpiarg" ;;
          -l*) MPI_LDFLAGS="$MPI_LDFLAGS $mpiarg" ;;
          -pthread) MPI_LDFLAGS="$MPI_LDFLAGS -lpthread" ;;
          -Wl*) MPI_CFLAGS="$MPI_CFLAGS $mpiarg" ;;
        esac
      done
    fi
  fi
  AC_MSG_RESULT([$MPI_CPPFLAGS $MPI_CFLAGS $MPI_LDFLAGS])
  AC_SUBST(MPI_TYPE)dnl
])


dnl This is AC_CHECK_SIZEOF but prepends LAL.

AC_DEFUN([LAL_CHECK_SIZEOF],
[changequote(<<, >>)dnl
dnl The name to #define.
define(<<LAL_TYPE_NAME>>, translit(lal_sizeof_$1, [a-z *], [A-Z_P]))dnl
dnl The cache variable name.
define(<<LAL_CV_NAME>>, translit(lal_cv_sizeof_$1, [ *], [_p]))dnl
changequote([, ])dnl
AC_MSG_CHECKING([size of $1])
AC_CACHE_VAL(LAL_CV_NAME,
[AC_TRY_RUN([#include <stdio.h>
main()
{
  FILE *f=fopen("conftestval", "w");
  if (!f) exit(1);
  fprintf(f, "%d\n", sizeof($1));
  exit(0);
}], LAL_CV_NAME=`cat conftestval`, LAL_CV_NAME=0, ifelse([$2], , , LAL_CV_NAME=$2))
])dnl
AC_MSG_RESULT([$LAL_CV_NAME])
AC_DEFINE_UNQUOTED(LAL_TYPE_NAME, $LAL_CV_NAME)
undefine([LAL_TYPE_NAME])dnl
undefine([LAL_CV_NAME])dnl
])
