dnl acinclude.m4

AC_DEFUN(LALAPPS_WITH_EXTRA_CPPFLAGS,
[AC_ARG_WITH(
	extra_cppflags, 
        [  --with-extra-cppflags=CPPFLAGS  additional C preprocessor flags],
	[ if test -n "${with_extra_cppflags}"
	  then
	    CPPFLAGS="$CPPFLAGS ${with_extra_cppflags}";
	  fi
	],)
])

AC_DEFUN(LALAPPS_WITH_EXTRA_CFLAGS,
[AC_ARG_WITH(
	extra_cflags, 
        [  --with-extra-cflags=CFLAGS  additional C compiler flags],
	[ if test -n "${with_extra_cflags}"
	  then
	    CFLAGS="$CFLAGS ${with_extra_cflags}";
	  fi
	],)
])

AC_DEFUN(LALAPPS_WITH_EXTRA_LDFLAGS,
[AC_ARG_WITH(
	extra_ldflags, 
        [  --with-extra-ldflags=LDFLAGS  additional linker flags],
	[ if test -n "${with_extra_ldflags}"
	  then
	    LDFLAGS="$LDFLAGS ${with_extra_ldflags}";
	  fi
	],)
])

AC_DEFUN(LALAPPS_WITH_EXTRA_LIBS,
[AC_ARG_WITH(
	extra_libs, 
        [  --with-extra-libs=LIBS  additional -l and -L linker flags],
	[ if test -n "${with_extra_libs}"
	  then
	    LIBS="$LIBS ${with_extra_libs}";
	  fi
	],)
])

AC_DEFUN(LALAPPS_WITH_CC,
[AC_ARG_WITH(
        cc, 
        [  --with-cc=CC            use the CC C compiler],
        [ if test -n "${with_cc}"
          then
            CC="${with_cc}";
          fi
        ],)
])

AC_DEFUN(LALAPPS_CHECK_LAL,
[AC_MSG_CHECKING([for -llal])
AC_CACHE_VAL(ac_cv_lib_lal,
[ac_save_LIBS="$LIBS"
LIBS="-llal $LIBS"
AC_TRY_LINK([int lalDebugLevel = 0; char LALMalloc();], LALMalloc();,
eval "ac_cv_lib_lal=yes", eval "ac_cv_lib_lal=no")
LIBS="$ac_save_LIBS"
])dnl
if eval "test \"`echo '$ac_cv_lib_lal'`\" = yes"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE(HAVE_LIBLAL, 1, [Define if you have the lal library (-llal).])
  LIBS="-llal $LIBS"
else
  AC_MSG_RESULT(no)
  AC_MSG_ERROR(could not find the LAL library)
fi
])







