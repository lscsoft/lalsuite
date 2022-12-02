# -*- mode: autoconf; -*-
# lalsuite_build.m4 - top level build macros
#
# serial 170

# restrict which LALSUITE_... patterns can appearing in output (./configure);
# useful for debugging problems with unexpanded LALSUITE_... Autoconf macros
m4_pattern_forbid([^_?LALSUITE_[A-Z_]+$])
m4_pattern_allow([^LALSUITE_(BUILD)$])

# list of user variables; see section 4.8.1 of the Autoconf manual
m4_define([uvar_list],[CPPFLAGS CFLAGS CXXFLAGS FCFLAGS FFLAGS LDFLAGS])
# prefix used to save/restore user variables in
m4_define([uvar_orig_prefix],[lalsuite_uvar_])
m4_define([uvar_prefix],uvar_orig_prefix)

AC_DEFUN([LALSUITE_ARG_VAR],[
  AC_ARG_VAR(LAL_BIN_PATH,[Location of LAL executables])
  AC_ARG_VAR(LAL_DATA_PATH,[Location of LAL data files])
  AC_ARG_VAR(LAL_OCTAVE_PATH,[Location of LAL octave files])
  AC_ARG_VAR(LAL_PYTHON_PATH,[Location of LAL python files])
  AC_ARG_VAR(LALSUITE_BUILD,[Set if part of lalsuite build])
])

m4_append([AC_INIT],[
  # just after AC_INIT:
  # save user-supplied values of user variables
  m4_foreach_w([uvar],uvar_list,[
    uvar_prefix[]uvar="${uvar}"
  ])
  m4_pushdef([uvar_prefix],uvar_prefix[]p_)
])

m4_rename([AC_OUTPUT],[lalsuite_AC_OUTPUT])
AC_DEFUN([AC_OUTPUT],[
  # just before AC_OUTPUT:
  # check for unbalanced LALSUITE_{PUSH,POP}_UVAR pairs
  m4_popdef([uvar_prefix])
  m4_if(uvar_prefix,uvar_orig_prefix,[],[
    m4_fatal([unbalanced LALSUITE_{PUSH,POP}_UVAR pairs])
  ])
  # prepend compiler configuration e.g. CFLAGS to AM_CFLAGS,
  # then restore original user-supplied values of user variables
  m4_foreach_w([uvar],uvar_list,[
    AC_SUBST(AM_[]uvar,"${AM_[]uvar} ${sys_[]uvar}")
    uvar="${uvar_prefix[]uvar}"
  ])
  # append extra values for user variables to be added after configuration
  m4_foreach_w([uvar],uvar_list,[
    AC_ARG_VAR([POSTCONFIG_]uvar,[Extra ]uvar[ to be added after configuration])
    uvar="${uvar} ${POSTCONFIG_[]uvar}"
  ])
  # call original AC_OUTPUT
  lalsuite_AC_OUTPUT
])

AC_DEFUN([LALSUITE_PUSH_UVARS],[
  # $0: save current values of user variables and LIBS
  m4_foreach_w([uvar],uvar_list[ LIBS],[
    uvar_prefix[]uvar="${uvar}"
    uvar_prefix[]uvar[]_lineno=$LINENO; _AS_ECHO_LOG([pushed uvar=${uvar}])
  ])
  m4_pushdef([uvar_prefix],uvar_prefix[]p_)
  # end $0
])

AC_DEFUN([LALSUITE_CLEAR_UVARS],[
  # $0: clear current values of user variables and LIBS
  m4_foreach_w([uvar],uvar_list[ LIBS],[
    AS_UNSET(uvar)
  ])
  # end $0
])

AC_DEFUN([LALSUITE_POP_UVARS],[
  # $0: restore previous values of user variables and LIBS
  m4_popdef([uvar_prefix])
  m4_foreach_w([uvar],uvar_list[ LIBS],[
    uvar="${uvar_prefix[]uvar}"
   _AS_ECHO_LOG([popped uvar from line ${uvar_prefix[]uvar[]_lineno}])
  ])
  # end $0
])

AC_DEFUN([LALSUITE_CHECK_PLATFORM],[
  # $0: detect the build platform and set conditionals
  host_darwin=no
  case $host_os in
    darwin*)
      host_darwin=yes
      ;;
  esac
  AM_CONDITIONAL([HOST_OS_DARWIN], [test x"${host_darwin}" = xyes])
  AM_COND_IF([HOST_OS_DARWIN],[
    LALSUITE_ADD_TESTS_ENV_CONFIG_VAR([HOST_OS_DARWIN],[true])
  ],[
    LALSUITE_ADD_TESTS_ENV_CONFIG_VAR([HOST_OS_DARWIN],[false])
  ])
])

AC_DEFUN([LALSUITE_CHECK_COMPILE_FLAGS],[
  # $0: check multiple compile flags
  LALSUITE_PUSH_UVARS
  LALSUITE_CLEAR_UVARS
  lalsuite_save_werror_flag=${ac_[]_AC_LANG_ABBREV[]_werror_flag}
  ac_[]_AC_LANG_ABBREV[]_werror_flag=yes
  lalsuite_passed_flags=""
  for flag in m4_normalize($1); do
    AS_VAR_PUSHDEF([CACHEVAR],[lalsuite_cv_check_compile_[]_AC_LANG_ABBREV[]_${flag}])
    AC_CACHE_CHECK([whether ]_AC_LANG[ compiler accepts ${flag}],CACHEVAR,[
      _AC_LANG_PREFIX[]FLAGS="${uvar_orig_prefix[]_AC_LANG_PREFIX[]FLAGS} ${lalsuite_passed_flags} ${flag}"
      AC_COMPILE_IFELSE([AC_LANG_SOURCE([])],[
        AS_VAR_SET(CACHEVAR,[yes])
      ],[
        AS_VAR_SET(CACHEVAR,[no])
      ])
    ])
    AS_IF([test x"AS_VAR_GET(CACHEVAR)" = xyes],[
      m4_default([$2], :)
      lalsuite_passed_flags="${lalsuite_passed_flags} ${flag}"
    ],[
      m4_default([$3], :)
    ])
    AS_VAR_POPDEF([CACHEVAR])
  done
  ac_[]_AC_LANG_ABBREV[]_werror_flag=${lalsuite_save_werror_flag}
  LALSUITE_POP_UVARS
  # end $0
])

AC_DEFUN([LALSUITE_CHECK_LINK_FLAGS],[
  # $0: check multiple link flags
  LALSUITE_PUSH_UVARS
  LALSUITE_CLEAR_UVARS
  lalsuite_save_werror_flag=${ac_[]_AC_LANG_ABBREV[]_werror_flag}
  ac_[]_AC_LANG_ABBREV[]_werror_flag=yes
  lalsuite_passed_flags=""
  for flag in m4_normalize($1); do
    AS_VAR_PUSHDEF([CACHEVAR],[lalsuite_cv_check_link_[]_AC_LANG_ABBREV[]_${flag}])
    AC_CACHE_CHECK([whether ]_AC_LANG[ linker accepts ${flag}],CACHEVAR,[
      LDFLAGS="${uvar_orig_prefix[]LDFLAGS} ${lalsuite_passed_flags} ${flag}"
      AC_LINK_IFELSE([AC_LANG_PROGRAM([])],[
        AS_VAR_SET(CACHEVAR,[yes])
      ],[
        AS_VAR_SET(CACHEVAR,[no])
      ])
    ])
    AS_IF([test x"AS_VAR_GET(CACHEVAR)" = xyes],[
      m4_default([$2],[:])
      lalsuite_passed_flags="${lalsuite_passed_flags} ${flag}"
    ],[
      m4_default([$3],[:])
    ])
    AS_VAR_POPDEF([CACHEVAR])
  done
  ac_[]_AC_LANG_ABBREV[]_werror_flag=${lalsuite_save_werror_flag}
  LALSUITE_POP_UVARS
  # end $0
])

AC_DEFUN([LALSUITE_ADD_FLAGS],[
  # $0: prepend/append flags to AM_CPPFLAGS/AM_$1FLAGS/AM_LDFLAGS/LIBS,
  # and update values of CPPFLAGS/$1FLAGS/LDFLAGS for Autoconf tests
  # - arg 1: prefix of the compiler flag variable, e.g. C for CFLAGS, CXX for CXXFLAGS
  # - arg 2: compiler flags
  # - arg 3: linker flags
  m4_ifval([$1],[m4_ifval([$2],[
    pre_AM_CPPFLAGS=
    post_AM_CPPFLAGS=
    pre_sys_CPPFLAGS=
    post_AM_$1FLAGS=
    for flag in $2; do
      # AM_CPPFLAGS gets unique -I, -D and -U flags
      # sys_CPPFLAGS gets unique system -I flags
      # AM_$1FLAGS gets everything else
      AS_CASE([${flag}],
        [-I/opt/*|-I/usr/*],[
          AS_CASE([" ${sys_CPPFLAGS} "],
            [*" ${flag} "*],[:],
            [pre_sys_CPPFLAGS="${pre_sys_CPPFLAGS} ${flag}"]
          )
        ],
        [-I*],[
          AS_CASE([" ${AM_CPPFLAGS} "],
            [*" ${flag} "*],[:],
            [pre_AM_CPPFLAGS="${pre_AM_CPPFLAGS} ${flag}"]
          )
        ],
        [-D*|-U*],[post_AM_CPPFLAGS="${post_AM_CPPFLAGS} ${flag}"],
        [post_AM_$1FLAGS="${post_AM_$1FLAGS} ${flag}"]
      )
    done
    AS_IF([test "x${pre_AM_CPPFLAGS}" != x],[
      AM_CPPFLAGS="${pre_AM_CPPFLAGS} ${AM_CPPFLAGS}"
      _AS_ECHO_LOG([prepended ${pre_AM_CPPFLAGS} to AM_CPPFLAGS])
    ])
    AS_IF([test "x${post_AM_CPPFLAGS}" != x],[
      AM_CPPFLAGS="${AM_CPPFLAGS} ${post_AM_CPPFLAGS}"
      _AS_ECHO_LOG([appended ${post_AM_CPPFLAGS} to AM_CPPFLAGS])
    ])
    AS_IF([test "x${pre_sys_CPPFLAGS}" != x],[
      sys_CPPFLAGS="${pre_sys_CPPFLAGS} ${sys_CPPFLAGS}"
      _AS_ECHO_LOG([prepended ${pre_sys_CPPFLAGS} to system AM_CPPFLAGS])
    ])
    AS_IF([test "x${post_AM_$1FLAGS}" != x],[
      AM_$1FLAGS="${AM_$1FLAGS} ${post_AM_$1FLAGS}"
      _AS_ECHO_LOG([appended ${post_AM_$1FLAGS} to AM_$1FLAGS])
    ])
    CPPFLAGS="${AM_CPPFLAGS} ${sys_CPPFLAGS} ${uvar_orig_prefix[]CPPFLAGS}"
    $1FLAGS="${AM_$1FLAGS} ${uvar_orig_prefix[]$1FLAGS}"
  ])])
  m4_ifval([$3],[
    pre_AM_LDFLAGS=
    post_AM_LDFLAGS=
    pre_sys_LDFLAGS=
    pre_LIBS=
    for flag in $3; do
      # LIBS gets -l flags and .la files
      # sys_LDFLAGS gets unique system -L flags
      # AM_LDFLAGS gets unique -L flags and everything else
      AS_CASE([${flag}],
        [-L/opt/*|-L/usr/*],[
          AS_CASE([" ${sys_LDFLAGS} "],
            [*" ${flag} "*],[:],
            [pre_sys_LDFLAGS="${pre_sys_LDFLAGS} ${flag}"]
          )
        ],
        [-L*],[
          AS_CASE([" ${AM_LDFLAGS} "],
            [*" ${flag} "*],[:],
            [
              flag_libdir=`expr "X${flag}" : "X-L\(.*\)"`
              pre_AM_LDFLAGS="${pre_AM_LDFLAGS} -L${flag_libdir} -R${flag_libdir}"
            ]
          )
        ],
        [-l*|*.la],[pre_LIBS="${pre_LIBS} ${flag}"],
        [post_AM_LDFLAGS="${post_AM_LDFLAGS} ${flag}"]
      )
    done
    AS_IF([test "x${pre_AM_LDFLAGS}" != x],[
      AM_LDFLAGS="${pre_AM_LDFLAGS} ${AM_LDFLAGS}"
      _AS_ECHO_LOG([prepended ${pre_AM_LDFLAGS} to AM_LDFLAGS])
    ])
    AS_IF([test "x${post_AM_LDFLAGS}" != x],[
      AM_LDFLAGS="${AM_LDFLAGS} ${post_AM_LDFLAGS}"
      _AS_ECHO_LOG([appended ${post_AM_LDFLAGS} to AM_LDFLAGS])
    ])
    AS_IF([test "x${pre_sys_LDFLAGS}" != x],[
      sys_LDFLAGS="${pre_sys_LDFLAGS} ${sys_LDFLAGS}"
      _AS_ECHO_LOG([prepended ${pre_sys_LDFLAGS} to system AM_LDFLAGS])
    ])
    AS_IF([test "x${pre_LIBS}" != x],[
      LIBS="${pre_LIBS} ${LIBS}"
      _AS_ECHO_LOG([prepended ${pre_LIBS} to LIBS])
    ])
    LDFLAGS="${AM_LDFLAGS} ${sys_LDFLAGS} ${uvar_orig_prefix[]LDFLAGS}"
  ])
  # end $0
])

AC_DEFUN([LALSUITE_ADD_PATH],[
  # $0: prepend path to $1, removing duplicates, first value taking precedence
  # - arg 1: name of path variable
  # - arg 2: path to prepend
  # - arg 3: whether to exclude /opt/... and /usr/... directories (default: yes)
  tokens=$2
  tokens=`echo ${tokens} ${$1} | sed 's/:/ /g'`
  $1=
  for token in ${tokens}; do
    AS_CASE([m4_default([$3],[yes]):${token}],
      [yes:/opt/*|yes:/usr/*],[:],
      AS_CASE([":${$1}:"],
        [*:${token}:*],[:],
        AS_IF([test "x${$1}" = x],[
          $1="${token}"
        ],[
          $1="${$1}:${token}"
        ])
      )
    )
  done
  _AS_ECHO_LOG([$1=${$1}])
  # end $0
])

AC_DEFUN([LALSUITE_VERSION_COMPARE],[
  # $0: compare versions using the given operator
  m4_case([$2],[<],,[<=],,[=],,[>=],,[>],,[m4_fatal([$0: invalid operator $2])])
  AS_VERSION_COMPARE([$1],[$3],[lalsuite_op='<'],[lalsuite_op='='],[lalsuite_op='>'])
  AS_CASE(['$2'],[*${lalsuite_op}*],[m4_default([$4],[:])],[m4_default([$5],[:])])
  # end $0
])

AC_DEFUN([LALSUITE_CHECK_GIT_REPO],[
  # $0: check for git
  AC_PATH_PROGS([GIT],[git],[false])
  # check whether building from a git repository
  have_git_repo=no
  AS_IF([test "x${GIT}" != xfalse],[
    AC_MSG_CHECKING([whether building from a git repository])
    # git log will print:
    # * the last log message, if the cwd is in a git repository
    # * nothing, if the cwd is not part of the git repo (e.g. ignored)
    # * an error msg to stderr if the cwd is not in a git repository
    git_log=`( cd "${srcdir}" && ${GIT} log --pretty=oneline -n 1 -- . ) 2>/dev/null`
    AS_IF([test "x${git_log}" != x],[have_git_repo=yes])
    AC_MSG_RESULT([${have_git_repo}])
  ])
  # conditional for git and building from a git repository
  AM_CONDITIONAL([HAVE_GIT_REPO],[test "x${have_git_repo}" = xyes])
  # end $0
])

AC_DEFUN([LALSUITE_VERSION_CONFIGURE_INFO],[
  # $0: define version/configure info
  m4_pushdef([uppercase],m4_translit(AC_PACKAGE_NAME, [a-z], [A-Z]))
  m4_pushdef([lowercase],m4_translit(AC_PACKAGE_NAME, [A-Z], [a-z]))
  m4_pushdef([withoutlal],m4_bpatsubst(AC_PACKAGE_NAME, [^LAL], []))
  version_major=`echo "$VERSION" | cut -d. -f1`
  version_minor=`echo "$VERSION" | cut -d. -f2`
  version_micro=`echo "$VERSION" | cut -d. -f3`
  version_devel=`echo "$VERSION" | cut -d. -f4-`
  test -z "$version_micro" && version_micro=0
  test -z "$version_devel" && version_devel=0
  reproducible_date="${srcdir}/gnuscripts/reproducible-date"
  configure_date=`${reproducible_date} +"%Y-%m-%dT%H:%M:%S%z"` >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD
  AS_IF([test "x${configure_date}" = x],[AC_MSG_ERROR([failed to run ${reproducible_date}])])
  AC_DEFINE_UNQUOTED(uppercase[_VERSION],["$VERSION"],AC_PACKAGE_NAME[ Version])
  AC_DEFINE_UNQUOTED(uppercase[_VERSION_MAJOR],[$version_major],AC_PACKAGE_NAME[ Version Major Number])
  AC_DEFINE_UNQUOTED(uppercase[_VERSION_MINOR],[$version_minor],AC_PACKAGE_NAME[ Version Minor Number])
  AC_DEFINE_UNQUOTED(uppercase[_VERSION_MICRO],[$version_micro],AC_PACKAGE_NAME[ Version Micro Number])
  AC_DEFINE_UNQUOTED(uppercase[_VERSION_DEVEL],[$version_devel],AC_PACKAGE_NAME[ Version Devel Number])
  AC_SUBST([ac_configure_args])
  AC_SUBST([configure_date])
  AC_SUBST([PACKAGE_NAME_UCASE],uppercase)
  AC_SUBST([PACKAGE_NAME_LCASE],lowercase)
  AC_SUBST([PACKAGE_NAME_NOLAL],withoutlal)
  m4_popdef([uppercase])
  m4_popdef([lowercase])
  m4_popdef([withoutlal])
  # end $0
])

AC_DEFUN([LALSUITE_ADD_TESTS_ENV_CONFIG_VAR],[
  AC_SUBST(TESTS_ENV_CONFIG_VARS)
  m4_if($2,[true],[
    AS_CASE([" ${TESTS_ENV_CONFIG_VARS} "],
      [*" export $1; $1=true; "*],[:],
      [*" export $1; $1=false; "*],[AC_MSG_ERROR([tried to add inconsistent value of $1=true (already $1=false) to TESTS_ENV_CONFIG_VARS])],
      [*],[
        TESTS_ENV_CONFIG_VARS="${TESTS_ENV_CONFIG_VARS} export $1; $1=true;"
        _AS_ECHO_LOG([added $1=true to TESTS_ENV_CONFIG_VARS])
      ]
    )
  ],[
    AS_CASE([" ${TESTS_ENV_CONFIG_VARS} "],
      [*" export $1; $1=true; "*],[AC_MSG_ERROR([tried to add inconsistent value of $1=false (already $1=true) to TESTS_ENV_CONFIG_VARS])],
      [*" export $1; $1=false; "*],[:],
      [*],[
        TESTS_ENV_CONFIG_VARS="${TESTS_ENV_CONFIG_VARS} export $1; $1=false;"
        _AS_ECHO_LOG([added $1=false to TESTS_ENV_CONFIG_VARS])
      ]
    )
  ])
])

AC_DEFUN([LALSUITE_ENABLE_STRICT_DEFS],[
  # $0: forbid certain definitions in LAL code
  AC_ARG_ENABLE(
    [strict-defs],
    AS_HELP_STRING(
      [--enable-strict-defs],
      [forbid certain definitions in LAL code [default=no]]
    ),[
      AS_CASE(["${enableval}"],
        [yes],[LALSUITE_ADD_FLAGS([C],[-DLAL_STRICT_DEFS_ENABLED])],
        [no],[:],
        [AC_MSG_ERROR([bad value for ${enableval} for --enable-strict-defs])]
      )
    ],[
    ]
  )
  # end $0
])

AC_DEFUN([LALSUITE_REQUIRE_CXX],[
  # $0: require a C++ compiler
  lalsuite_require_cxx=true
  # end $0
])

AC_DEFUN([LALSUITE_REQUIRE_F77],[
  # $0: require an F77 compiler
  lalsuite_require_f77=true
  # end $0
])

# because we want to conditionally decide whether to check for
# C++/Fortran compilers only at ./configure run time, we must erase the
# following macros; in Autoconf 2.64 and later, they AC_REQUIRE the
# C++/Fortran AC_PROG_... macros, which forces the C++/Fortran compilers
# to always be checked for, which prevents us from instead conditionally
# deciding that at ./configure run time
m4_foreach([lang],[[C++],[Fortran 77],[Fortran]],[
  m4_defun([AC_LANG_COMPILER(]lang[)],[])
  m4_defun([AC_LANG_PREPROC(]lang[)],[])
])

AC_DEFUN([LALSUITE_PROG_COMPILERS],[
  # $0: check for C/C++/Fortran compilers
  AC_REQUIRE([AC_PROG_CC])
  AC_REQUIRE([AC_PROG_CPP])

  # check for C99 compiler
  # NOTE: AC_PROG_CC_C99 is deprecated in autoconf-2.70, so avoid calling it;
  #       on systems that have 2.70 c99 or greater will be the compiler
  #       default anyway, so we don't _need_ to use AC_PROG_CC_C99
  m4_version_prereq([2.70], [], [AC_REQUIRE([AC_PROG_CC_C99])])

  # set default CFLAGS
  CFLAGS=
  AS_IF([test "x${GCC}" = xyes],[
    CFLAGS="${CFLAGS} -O2"
  ])
  AS_IF([test "x${ac_cv_prog_cc_g}" = xyes],[
    CFLAGS="${CFLAGS} -g"
  ])

  # only include ISO C99 definitions from system headers
  CFLAGS="${CFLAGS} -D_ISOC99_SOURCE"

  # check for C++ compiler, if needed
  AS_IF([test "${lalsuite_require_cxx}" = true],[
    AC_PROG_CXX
    AC_PROG_CXXCPP

    # set default CXXFLAGS
    CXXFLAGS=
    AS_IF([test "x${GXX}" = xyes],[
      CXXFLAGS="${CXXFLAGS} -O2"
    ])
    AS_IF([test "x${ac_cv_prog_cxx_g}" = xyes],[
      CXXFLAGS="${CXXFLAGS} -g"
    ])

    # define C99 constant, limit, and format macros for C++ sources
    CXXFLAGS="${CXXFLAGS} -D__STDC_CONSTANT_MACROS -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS"

  ],[
    CXX=
    CXXCPP=
    CXXFLAGS=
    AM_CONDITIONAL([am__fastdepCXX],[false])
  ])

  # check for F77 compiler, if needed
  AS_IF([test "${lalsuite_require_f77}" = true],[
    AC_PROG_F77

    # set default FFLAGS
    FFLAGS=
    AS_IF([test "x${ac_cv_f77_compiler_gnu}" = xyes],[
      FFLAGS="${FFLAGS} -O2"
    ])
    AS_IF([test "x${ac_cv_prog_f77_g}" = xyes],[
      FFLAGS="${FFLAGS} -g"
    ])

  ],[
    F77=
    FFLAGS=
  ])

  # add flags
  LALSUITE_ADD_FLAGS([C],[${CFLAGS}])
  LALSUITE_ADD_FLAGS([CXX],[${CXXFLAGS}])
  LALSUITE_ADD_FLAGS([F],[${FFLAGS}])

  # end $0
])

AC_DEFUN([LALSUITE_REQUIRE_PYTHON],[
  # $0: require Python version $1 or later
  AS_IF([test "x${lalsuite_require_pyvers}" = x],[
    lalsuite_require_pyvers="$1"
  ],[
    LALSUITE_VERSION_COMPARE([$1],[>],[${lalsuite_require_pyvers}],[
      lalsuite_require_pyvers="$1"
    ])
  ])
  # end $0
])

AC_DEFUN([LALSUITE_CHECK_PYTHON],[
  # $0: check for Python
  AC_ARG_ENABLE(
    [python],
    AS_HELP_STRING(
      [--enable-python],
      [build Python programs and modules]
    ),[
      AS_CASE(["${enableval}"],
        [yes],[python=true],
        [no],[python=false],
        [AC_MSG_ERROR([invalid value "${enableval}" for --enable-python])]
      )
    ],[
      python=
    ]
  )
  lalsuite_pyvers="$1"
  AS_IF([test "x${lalsuite_require_pyvers}" != x],[
    LALSUITE_VERSION_COMPARE([${lalsuite_require_pyvers}],[>],[${lalsuite_pyvers}],[
      lalsuite_pyvers="${lalsuite_require_pyvers}"
    ])
  ])
  AS_IF([test "x${python}" != xfalse || test "x${lalsuite_require_pyvers}" != x],[
    AM_PATH_PYTHON([${lalsuite_pyvers}],[
      AC_SUBST([python_prefix], [`${PYTHON} -c 'import sys; print(sys.prefix)' 2>/dev/null`])
      AC_SUBST([python_exec_prefix], [`${PYTHON} -c 'import sys; print(sys.exec_prefix)' 2>/dev/null`])
    ],[
      AS_IF([test "x${python}" = xtrue || test "x${lalsuite_require_pyvers}" != x],[
        AC_MSG_ERROR([Python version ${lalsuite_pyvers} or later is required])
      ],[
        python=false
      ])
    ])
  ])
  AS_IF([test "x${python}" = xfalse && test "x${lalsuite_require_pyvers}" = x],[
    AC_SUBST([PYTHON],["${SHELL} -c 'echo ERROR: Python is missing >&2; exit 1' --"])
  ])
  AM_CONDITIONAL([HAVE_PYTHON],[test "x${python}" != xfalse])
  AM_COND_IF([HAVE_PYTHON],[
    PYTHON_ENABLE_VAL=ENABLED
    LALSUITE_ADD_TESTS_ENV_CONFIG_VAR([HAVE_PYTHON],[true])
  ],[
    PYTHON_ENABLE_VAL=DISABLED
    LALSUITE_ADD_TESTS_ENV_CONFIG_VAR([HAVE_PYTHON],[false])
  ])
  # end $0
])

AC_DEFUN([LALSUITE_USE_LIBTOOL],[
  # $0: Generate a libtool script for use in configure tests. Arguments
  # are added to link command in variable ${lalsuite_libtool_flags}
  AC_REQUIRE([LT_INIT])
  LT_OUTPUT
  m4_append([AC_LANG(C)],[
    ac_link="./libtool --mode=link --tag=CC ${ac_link} ${lalsuite_libtool_flags}"
  ])
  AC_PROVIDE_IFELSE([AC_PROG_CXX],[
    m4_append([AC_LANG(C++)],[
      ac_link="./libtool --mode=link --tag=CXX ${ac_link} ${lalsuite_libtool_flags}"
    ])
  ])
  AC_LANG(_AC_LANG)
  LALSUITE_ADD_FLAGS([],[],[${lalsuite_libtool_flags}])
  # end $0
])

AC_DEFUN([LALSUITE_DISTCHECK_CONFIGURE_FLAGS],[
  # $0: store configure flags for 'make distcheck'
  DISTCHECK_CONFIGURE_FLAGS=
  eval set x "${ac_configure_args}"
  shift
  for arg
  do
    AS_CASE(["${arg}"],
      [--enable-*|--disable-*],[
        # save any --enable/--disable arguments
        DISTCHECK_CONFIGURE_FLAGS="${DISTCHECK_CONFIGURE_FLAGS} '${arg}'"
      ],
      [--with-*|--without-*],[
        # save any --with/--without arguments
        DISTCHECK_CONFIGURE_FLAGS="${DISTCHECK_CONFIGURE_FLAGS} '${arg}'"
      ],
      [--*],[
        # skip all other ./configure arguments
        :
      ],
      [DISTCHECK_CONFIGURE_FLAGS=*],[
        # append value of DISTCHECK_CONFIGURE_FLAGS
        DISTCHECK_CONFIGURE_FLAGS="${DISTCHECK_CONFIGURE_FLAGS} '`expr "X${arg}" : "XDISTCHECK_CONFIGURE_FLAGS=\(.*\)"`'"
      ],
      [*=*],[
        # save any environment variables given to ./configure
        DISTCHECK_CONFIGURE_FLAGS="${DISTCHECK_CONFIGURE_FLAGS} '${arg}'"
      ]
    )
  done
  AC_SUBST(DISTCHECK_CONFIGURE_FLAGS)
  # end $0
])

AC_DEFUN([LALSUITE_ENABLE_MODULE],[
  # $0: enable/disable module $1
  m4_pushdef([lowercase],m4_translit([[$1]], [A-Z], [a-z]))
  m4_pushdef([uppercase],m4_translit([[$1]], [a-z], [A-Z]))
  AM_CONDITIONAL(uppercase,[test "x${lowercase}" = xtrue])
  AS_IF([test "${lowercase}" = "true"],[
    uppercase[]_ENABLE_VAL=ENABLED
    LALSUITE_ADD_TESTS_ENV_CONFIG_VAR(uppercase[]_ENABLED,[true])
  ],[
    uppercase[]_ENABLE_VAL=DISABLED
    LALSUITE_ADD_TESTS_ENV_CONFIG_VAR(uppercase[]_ENABLED,[false])
  ])
  _AS_ECHO_LOG([module $1 is ${]uppercase[_ENABLE_VAL}])
  m4_popdef([lowercase])
  m4_popdef([uppercase])
  # end $0
])

AC_DEFUN([LALSUITE_CHECK_LIB],[
  # $0: check for LAL library
  # - arg 1: name of LAL library
  # - arg 2: minimum version required
  # - arg 3: library function to check for
  # - arg 4: library header to check for
  m4_pushdef([lowercase],m4_translit([[$1]], [A-Z], [a-z]))
  m4_pushdef([uppercase],m4_translit([[$1]], [a-z], [A-Z]))

  # substitute required library version in pkg-config files
  AC_SUBST(uppercase[]_VERSION,[$2])

  # set up pkg-config environment
  AS_UNSET([PKG_CONFIG_ALLOW_SYSTEM_CFLAGS])
  AS_UNSET([PKG_CONFIG_ALLOW_SYSTEM_LIBS])

  # prepend to CFLAGS, CPPFLAGS, LDFLAGS, LIBS, LAL_BIN_PATH, LAL_DATA_PATH, LAL_OCTAVE_PATH, LAL_PYTHON_PATH, LAL_HTMLDIR
  PKG_CHECK_MODULES(uppercase, [lowercase >= $2], [lowercase="true"], [lowercase="false"])
  PKG_CHECK_VAR(uppercase[]_BIN_PATH, [lowercase >= $2], uppercase[]_BIN_PATH,,)
  PKG_CHECK_VAR(uppercase[]_DATA_PATH, [lowercase >= $2], uppercase[]_DATA_PATH,,)
  PKG_CHECK_VAR(uppercase[]_OCTAVE_PATH, [lowercase >= $2], uppercase[]_OCTAVE_PATH,,)
  PKG_CHECK_VAR(uppercase[]_PYTHON_PATH, [lowercase >= $2], uppercase[]_PYTHON_PATH,,)
  PKG_CHECK_VAR(uppercase[]_HTMLDIR, [lowercase >= $2], htmldir,,)
  if test "$lowercase" = "true"; then
    LALSUITE_ADD_FLAGS([C],$[]uppercase[]_CFLAGS,$[]uppercase[]_LIBS)
    LALSUITE_ADD_PATH(LAL_BIN_PATH,"$[]uppercase[]_BIN_PATH")
    LALSUITE_ADD_PATH(LAL_DATA_PATH,"$[]uppercase[]_DATA_PATH")
    LALSUITE_ADD_PATH(LAL_OCTAVE_PATH,"$[]uppercase[]_OCTAVE_PATH")
    LALSUITE_ADD_PATH(LAL_PYTHON_PATH,"$[]uppercase[]_PYTHON_PATH")
  fi

  # add system include flags to LAL_SYSTEM_INCLUDES
  if test -n "$PKG_CONFIG" -a "$LALSUITE_BUILD" != "true"; then
    # use pkg-config to get system paths
    PKG_CONFIG_ALLOW_SYSTEM_CFLAGS=1
    export PKG_CONFIG_ALLOW_SYSTEM_CFLAGS
    for flag in `$PKG_CONFIG --cflags-only-I "lowercase >= $2"`; do
      AS_CASE([" $CPPFLAGS $LAL_SYSTEM_INCLUDES "],
        [*" ${flag} "*],[:],
        [LAL_SYSTEM_INCLUDES="$LAL_SYSTEM_INCLUDES $flag"]
      )
    done
    AS_UNSET([PKG_CONFIG_ALLOW_SYSTEM_CFLAGS])
  else
    # use standard include paths
    save_IFS="$IFS"
    IFS=:
    for flag in $C_INCLUDE_PATH $CPLUS_INCLUDE_PATH ; do
      test -n "$flag" && flag="-I$flag"
      AS_CASE([" $CPPFLAGS $LAL_SYSTEM_INCLUDES "],
        [*" ${flag} "*],[:],
        [LAL_SYSTEM_INCLUDES="$LAL_SYSTEM_INCLUDES $flag"]
      )
    done
    IFS="$save_IFS"
  fi
  AC_SUBST([LAL_SYSTEM_INCLUDES])

  if test "$LALSUITE_BUILD" = "true"; then
    if test "$lowercase" = "false"; then
      # should never get here: bug in build system
      AC_MSG_ERROR([could not find the $1 library])
    fi
  else
    AC_CHECK_LIB(lowercase,[$3],,[AC_MSG_ERROR([could not find the $1 library])])
    AC_CHECK_HEADERS([$4],,[AC_MSG_ERROR([could not find the $4 header])])
    m4_if(lowercase,[lalsupport],[],[
      LALSUITE_HEADER_LIBRARY_MISMATCH_CHECK([$1])
    ])
  fi
  AC_DEFINE([HAVE_LIB]uppercase,[1],[Define to 1 if you have the $1 library])
  # add $1 to list of LALSuite libraries
  lalsuite_libs="${lalsuite_libs} lowercase"
  lowercase="true"
  LALSUITE_ENABLE_MODULE($1)
  m4_popdef([lowercase])
  m4_popdef([uppercase])
  # end $0
])

AC_DEFUN([LALSUITE_CHECK_OPT_LIB],[
  # $0: check for optional LAL library
  # - arg 1: name of LAL library
  # - arg 2: minimum version required
  # - arg 3: library function to check for
  # - arg 4: library header to check for
  m4_pushdef([lowercase],m4_translit([[$1]], [A-Z], [a-z]))
  m4_pushdef([uppercase],m4_translit([[$1]], [a-z], [A-Z]))

  # optional check for $1
  AS_IF([test "${lowercase}" = "true"],[
    LALSUITE_CHECK_LIB($1,$2,$3,$4)
  ])

  # enable/disable $1
  LALSUITE_ENABLE_MODULE($1)

  m4_popdef([lowercase])
  m4_popdef([uppercase])
  # end $0
])

AC_DEFUN([LALSUITE_HEADER_LIBRARY_MISMATCH_CHECK],[
  # $0: check for version mismatch between library $1 and its headers
  m4_pushdef([uppercase],m4_translit([[$1]], [a-z], [A-Z]))
  m4_pushdef([withoutlal],m4_bpatsubst([$1], [^LAL], []))
  AC_MSG_CHECKING([whether $1 headers match the library])
  AC_LINK_IFELSE([
    AC_LANG_SOURCE([[
#include <lal/$1VCSInfoHeader.h>
int main(void) {
#if ]uppercase[_VERSION_DEVEL != 0
  ]uppercase[_VCS_LINK_CHECK();
#endif
  return 0;
}
    ]])
  ],[
    AC_MSG_RESULT([yes])
  ],[
    AC_MSG_RESULT([no])
    AC_LINK_IFELSE([
      AC_LANG_SOURCE([[
#include <stdio.h>
#include <lal/$1VCSInfoHeader.h>
int main(void) {
  printf("$1 headers: %s %s %s\n", ]uppercase[_VERSION, ]uppercase[_VCS_ID, ]uppercase[_VCS_STATUS);
  printf("$1 library: %s %s %s\n", lal]withoutlal[VCSInfo.version, lal]withoutlal[VCSInfo.vcsId, lal]withoutlal[VCSInfo.vcsStatus);
  return 0;
}
      ]])
    ],[
      AS_IF([test ${cross_compiling} != yes],[
        ./conftest${EXEEXT} | ${SED} -e "s/^/${as_me}:${as_lineno-$LINENO}: /" >&AS_MESSAGE_LOG_FD
      ])
    ],[
      _AS_ECHO_LOG([could not determine further details of $1 header-library mismatch (compile/link error)])
    ],[
      _AS_ECHO_LOG([could not determine further details of $1 header-library mismatch (cross compiling)])
    ])
    AC_MSG_ERROR([Your $1 headers do not match your library. Check config.log for details.])
  ])
  m4_popdef([uppercase])
  m4_popdef([withoutlal])
  # end $0
])

AC_DEFUN([LALSUITE_CHECK_LIBRARY_FOR_SUPPORT],[
  # $0: check if library $1 supports $2
  m4_pushdef([uppercase],m4_translit([[$1]], [a-z], [A-Z]))
  LALSUITE_PUSH_UVARS
  # because library may be uninstalled, remove everything but -I flags
  save_CPPFLAGS="${CPPFLAGS}"
  LALSUITE_CLEAR_UVARS
  CPPFLAGS="${save_CPPFLAGS}"
  AS_IF([test "x${LALSUITE_BUILD}" = xtrue],[
    for arg in ${lalsuite_libs}; do
      AS_CASE([${arg}],
        [lalsupport],[:],[
          CPPFLAGS="-I${ac_pwd}/../${arg}/src ${CPPFLAGS}"
        ]
      )
    done
    CPPFLAGS="-DLALSUITE_BUILD ${CPPFLAGS}"
  ])
  AC_MSG_CHECKING([whether $1 has been compiled with $2 support])
  AC_COMPILE_IFELSE([
    AC_LANG_SOURCE([[
#ifdef LALSUITE_BUILD
#include "$1Config.h"
#else
#include <lal/$1Config.h>
#endif
#ifndef ]uppercase[_$2_ENABLED
#error ]uppercase[_$2_ENABLED is not defined
#endif
    ]])
  ],[
    AC_MSG_RESULT([yes])
    m4_default([$3],[:])
  ],[
    AC_MSG_RESULT([no])
    m4_default([$4],[:])
  ])
  LALSUITE_POP_UVARS
  m4_popdef([uppercase])
  # end $0
])

AC_DEFUN([LALSUITE_ENABLE_NIGHTLY],
[
  BASE_VERSION="${VERSION}"
  AC_ARG_ENABLE(
    [nightly],
    AS_HELP_STRING([--enable-nightly],[nightly build [default=no]]),
    [ case "${enableval}" in
        yes) NIGHTLY_VERSION=dev`date -u +"%Y%m%d"`
             VERSION="${BASE_VERSION}-${NIGHTLY_VERSION}" ;;
        no) NIGHTLY_VERSION="";;
        *) NIGHTLY_VERSION="${enableval}"
           VERSION="${BASE_VERSION}-${NIGHTLY_VERSION}" ;;
        esac ],
    [ NIGHTLY_VERSION="" ] )
    AC_SUBST([VERSION])
    AC_SUBST([BASE_VERSION])
    AC_SUBST([NIGHTLY_VERSION])
])

AC_DEFUN([LALSUITE_ENABLE_ALL_LAL],
[AC_ARG_ENABLE(
  [all_lal],
  AS_HELP_STRING([--enable-all-lal],[enable/disable compilation of all LAL libraries]),
  [ case "${enableval}" in
      yes) all_lal=true;;
      no) all_lal=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-all-lal) ;;
    esac
  ], [ all_lal= ] )
])

AC_DEFUN([LALSUITE_ENABLE_LALFRAME],
[AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
AC_ARG_ENABLE(
  [lalframe],
  AS_HELP_STRING([--enable-lalframe],[compile code that requires lalframe library [default=yes]]),
  [ case "${enableval}" in
      yes) lalframe=true;;
      no) lalframe=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalframe) ;;
    esac
  ], [ lalframe=${all_lal:-true} ] )
])

AC_DEFUN([LALSUITE_ENABLE_LALMETAIO],
[AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
AC_ARG_ENABLE(
  [lalmetaio],
  AS_HELP_STRING([--enable-lalmetaio],[compile code that requires lalmetaio library [default=yes]]),
  [ case "${enableval}" in
      yes) lalmetaio=true;;
      no) lalmetaio=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalmetaio) ;;
    esac
  ], [ lalmetaio=${all_lal:-true} ] )
])

AC_DEFUN([LALSUITE_ENABLE_LALSIMULATION],
[AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
AC_ARG_ENABLE(
  [lalsimulation],
  AS_HELP_STRING([--enable-lalsimulation],[compile code that requires lalsimulation library [default=yes]]),
  [ case "${enableval}" in
      yes) lalsimulation=true;;
      no) lalsimulation=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalsimulation) ;;
    esac
  ], [ lalsimulation=${all_lal:-true} ] )
])

AC_DEFUN([LALSUITE_ENABLE_LALBURST],
[AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
AC_ARG_ENABLE(
  [lalburst],
  AS_HELP_STRING([--enable-lalburst],[compile code that requires lalburst library [default=yes]]),
  [ case "${enableval}" in
      yes) lalburst=true;;
      no) lalburst=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalburst) ;;
    esac
  ], [ lalburst=${all_lal:-true} ] )
])

AC_DEFUN([LALSUITE_ENABLE_LALINSPIRAL],
[AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
AC_ARG_ENABLE(
  [lalinspiral],
  AS_HELP_STRING([--enable-lalinspiral],[compile code that requires lalinspiral library [default=yes]]),
  [ case "${enableval}" in
      yes) lalinspiral=true;;
      no) lalinspiral=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalinspiral) ;;
    esac
  ], [ lalinspiral=${all_lal:-true} ] )
])

AC_DEFUN([LALSUITE_ENABLE_LALINFERENCE],
[AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
AC_ARG_ENABLE(
  [lalinference],
  AS_HELP_STRING([--enable-lalinference],[compile code that requires lalinference library [default=yes]]),
  [ case "${enableval}" in
      yes) lalinference=true;;
      no) lalinference=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalinference) ;;
    esac
  ], [ lalinference=${all_lal:-true} ] )
])

AC_DEFUN([LALSUITE_ENABLE_LALPULSAR],
[AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
AC_ARG_ENABLE(
  [lalpulsar],
  AS_HELP_STRING([--enable-lalpulsar],[compile code that requires lalpulsar library [default=yes]]),
  [ case "${enableval}" in
      yes) lalpulsar=true;;
      no) lalpulsar=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalpulsar) ;;
    esac
  ], [ lalpulsar=${all_lal:-true} ] )
])

AC_DEFUN([LALSUITE_ENABLE_LALAPPS],[
  AC_REQUIRE([LALSUITE_ENABLE_ALL_LAL])
  AC_ARG_ENABLE(
    [lalapps],
    AS_HELP_STRING([--enable-lalapps],[compile lalapps [default=yes]]),
    [
      case "${enableval}" in
        yes) lalapps=true ;;
        no) lalapps=false ;;
        *) AC_MSG_ERROR(bad value ${enableval} for --enable-lalapps) ;;
      esac
    ],[
      lalapps=${all_lal:-true}
    ]
  )
])

AC_DEFUN([LALSUITE_WITH_CUDA],[
  AC_ARG_WITH(
    [cuda],
    AS_HELP_STRING([--with-cuda=PATH],[specify location of CUDA [/opt/cuda]]),
    [
      AS_CASE([${with_cuda}],
        [no],[cuda=false],
        [yes],[cuda=true; cuda_path=/opt/cuda],
        [cuda=true; cuda_path=${with_cuda}]
      )
    ],[
      cuda=false
    ]
  )
  AC_ARG_WITH(
    [nvcc_cflags],
    AS_HELP_STRING([--with-nvcc-cflags=NVCC_CFLAGS],[NVCC compiler flags]),
    [
      NVCC_CFLAGS="$NVCC_CFLAGS ${with_nvcc_cflags}"
    ]
  )
  AS_IF([test "${cuda}" = true],[
    LALSUITE_REQUIRE_CXX
  ])
])

AC_DEFUN([LALSUITE_USE_CUDA],[
  AS_IF([test "${cuda}" = true],[
    AC_MSG_NOTICE([Using ${with_cuda} as CUDA path])
    AS_CASE([$build_os],
      [linux*],[
        AS_IF([test "x$build_cpu" = "xx86_64"],[
          cuda_libdir=lib64
        ],[
          cuda_libdir=lib
        ])
      ],
      [cuda_libdir=lib]
    )
    CUDA_LIBS="-L${cuda_path}/${cuda_libdir} -lcufft -lcudart"
    CUDA_CFLAGS="-I${with_cuda}/include"
    LALSUITE_ADD_FLAGS([C],[${CUDA_CFLAGS}],[${CUDA_LIBS}])
    AC_SUBST(CUDA_LIBS)
    AC_SUBST(CUDA_CFLAGS)
    AC_PATH_PROGS(NVCC,[nvcc],[],[${cuda_path}/bin:${PATH}])
    AS_IF([test "x${NVCC}" = x],[
      AC_MSG_ERROR([could not find 'nvcc' in path])
    ])
    AC_SUBST(NVCC_CFLAGS)
  ])
  LALSUITE_ENABLE_MODULE([CUDA])
])


AC_DEFUN([LALSUITE_CHECK_GSL_VERSION],[
  # $0: check for GSL version
  lal_min_gsl_version=m4_normalize([$1])
  AC_MSG_CHECKING(for GSL version >= $lal_min_gsl_version)
  AC_RUN_IFELSE(
    [AC_LANG_SOURCE([[
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
]])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_ERROR([could not find required version of GSL])],
    [AC_MSG_WARN([cross compiling; assumed OK...])])
  # end $0
])

AC_DEFUN([LALSUITE_ENABLE_FAST_GSL],[
  # $0: enable/disable fast/inline GSL code
  AC_ARG_ENABLE(
    [fast_gsl],
    AS_HELP_STRING([--enable-fast-gsl],[enable fast/inline GSL code [default=no]]),
    AS_CASE(["${enableval}"],
      [yes],[
        AC_DEFINE([HAVE_INLINE],[1],[Define to 1 to use inline code])
        AC_DEFINE([GSL_C99_INLINE],[1],[Define to 1 to use GSL C99 inline code])
        AC_DEFINE([GSL_RANGE_CHECK_OFF],[1],[Define to 1 to turn GSL range checking off])
      ],
      [no],[:],
      AC_MSG_ERROR([bad value ${enableval} for --enable-fast-gsl])
    )
  )
  # end $0
])

AC_DEFUN([LALSUITE_ENABLE_OSX_VERSION_CHECK],
[AC_ARG_ENABLE(
  [osx_version_check],
  AS_HELP_STRING([--enable-osx-version-check],[disable OS X version check [default=yes]]),
  AS_CASE(["${enableval}"],
    [yes],[osx_version_check=true],
    [no],[osx_version_check=false],
    AC_MSG_ERROR([bad value ${enableval} for --enable-osx-version-check])
  ),[osx_version_check=true])
])

AC_DEFUN([LALSUITE_OSX_VERSION_CHECK],[
LALSUITE_ENABLE_OSX_VERSION_CHECK
AS_IF([test "x${osx_version_check}" = "xtrue"],[
  AS_IF([test "x$build_vendor" = "xapple"],[
    AC_CHECK_PROGS([SW_VERS],[sw_vers])
    AS_IF([test "x$SW_VERS" != "x"],[
      AC_MSG_CHECKING([Mac OS X version])
      MACOSX_VERSION=`$SW_VERS -productVersion`
      AC_MSG_RESULT([$MACOSX_VERSION])])
    AS_CASE(["$MACOSX_VERSION"],
      [10.0*|10.1|10.1.*|10.2*|10.3*],AC_MSG_ERROR([This version of Mac OS X is not supported]),
      [10.4*|10.5*|10.6*|10.7*|10.8*|10.9*|10.10*|10.11*|10.12*|10.13*|10.14*],,
      AC_MSG_WARN([Unknown Mac OS X version]))
])])])

AC_DEFUN([LALSUITE_CHECK_CUDA],
[AC_MSG_CHECKING([whether LAL has been compiled with CUDA support])
AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <lal/LALConfig.h>
#ifdef LAL_CUDA_ENABLED
int main( void ) { return 0; }
#else
int main( void ) { return 1; }
#endif
]])],
AC_MSG_RESULT([yes])
[cuda=true],
AC_MSG_RESULT([no])
[cuda=false],
AC_MSG_RESULT([unknown])
[cuda=false])
])

AC_DEFUN([LALSUITE_ENABLE_DOXYGEN],[
  # $0: enable Doxygen documentation
  AC_ARG_ENABLE(
    [doxygen],
    AS_HELP_STRING(
      [--enable-doxygen],
      [generate Doxygen documentation]
    ),[
      AS_CASE(["${enableval}"],
        [yes],[doxygen=true],
        [no],[doxygen=false],
        [AC_MSG_ERROR([invalid value "${enableval}" for --enable-doxygen])]
      )
    ],[
      doxygen=false
    ]
  )
  # end $0
])

AC_DEFUN([LALSUITE_USE_DOXYGEN],[
  # $0: configure Doxygen documentation
  AC_REQUIRE([LALSUITE_CHECK_GIT_REPO])

  # add configuration option for MathJax installation
  AC_ARG_WITH(
    [mathjax],
    AS_HELP_STRING(
      [--with-mathjax=PATH],
      [use MathJax installation at PATH [default: use CDN]]
    ),[
      DOXYGEN_MATHJAXDIR="${withval}"
    ],[
      DOXYGEN_MATHJAXDIR=
    ]
  )

  # enable Doxygen module
  LALSUITE_ENABLE_MODULE([DOXYGEN])

  AS_IF([test "x${doxygen}" = xtrue],[

    # configure Doxygen filter script
    AC_CONFIG_FILES([doxygen/filter],[chmod +x doxygen/filter])

    # Python is required to run some scripts
    LALSUITE_REQUIRE_PYTHON([3.5])

    # Perl and BibTeX are required to build the references
    AC_PATH_PROG([PERL],[perl],[],[])
    AS_IF([test "x${PERL}" = x],[
      AC_MSG_ERROR([could not find 'perl' in PATH (required for Doxygen references)])
    ])
    AC_PATH_PROG([BIBTEX],[bibtex],[],[])
    AS_IF([test "x${BIBTEX}" = x],[
      AC_MSG_ERROR([could not find 'bibtex' in PATH (required for Doxygen references)])
    ])

    # check for Doxygen
    AC_PATH_PROG([DOXYGEN],[doxygen],[],[])
    AS_IF([test "x${DOXYGEN}" = x],[
      AC_MSG_ERROR([could not find 'doxygen' in PATH])
    ])
    doxygen_min_version=1.8.1.2   # minimum required Doxygen version
    AC_MSG_CHECKING([if ${DOXYGEN} version is at least ${doxygen_min_version}])
    doxygen_version=`${DOXYGEN} --version 2>/dev/null`
    LALSUITE_VERSION_COMPARE([${doxygen_version}],[<],[${doxygen_min_version}],[
      AC_MSG_RESULT([no (${doxygen_version})])
      AC_MSG_ERROR([Doxygen version ${doxygen_min_version} or later is required])
    ])
    AC_MSG_RESULT([yes (${doxygen_version})])

    # ignore some Doxygen warnings due to Doxygen bugs
    AC_SUBST([DOXYGEN_WARNING_REGEX],["-e '/^\$$/d'"])
    LALSUITE_VERSION_COMPARE([${doxygen_version}],[<],[1.8.9.1],[
      # https://bugzilla.gnome.org/show_bug.cgi?id=742151
      DOXYGEN_WARNING_REGEX=["${DOXYGEN_WARNING_REGEX} -e '/^citelist/d'"]
    ])
    LALSUITE_VERSION_COMPARE([1.8.8],[<=],[${doxygen_version}],[
      LALSUITE_VERSION_COMPARE([${doxygen_version}],[<],[1.8.11],[
        # https://bugzilla.gnome.org/show_bug.cgi?id=743604
        DOXYGEN_WARNING_REGEX=["${DOXYGEN_WARNING_REGEX} -e '/warning: Duplicate anchor/d'"]
        DOXYGEN_WARNING_REGEX=["${DOXYGEN_WARNING_REGEX} -e '/warning: multiple use of section label/d'"]
      ])
    ])
    LALSUITE_VERSION_COMPARE([1.8.8],[<=],[${doxygen_version}],[
      LALSUITE_VERSION_COMPARE([${doxygen_version}],[<=],[1.8.9.1],[
        # https://bugzilla.gnome.org/show_bug.cgi?id=743605
        DOXYGEN_WARNING_REGEX=["${DOXYGEN_WARNING_REGEX} -e '/warning: explicit link request/d'"]
      ])
    ])

    # build some substituted variables from list of configured LAL libraries
    AC_SUBST([DOXYGEN_ENABLED_SECTIONS])
    AC_SUBST([DOXYGEN_TAGFILES],[])
    AC_SUBST([DOXYGEN_INSTALL_DIRMAP],[])
    for arg in ${lalsuite_libs}; do
      AS_CASE([${arg}],
        [lalsupport],[:],[
          ucarg=[`echo ${arg} | tr '[a-z]' '[A-Z]'`]
          DOXYGEN_ENABLED_SECTIONS="${DOXYGEN_ENABLED_SECTIONS} ${ucarg}"
          arg_htmldir_name="${ucarg}_HTMLDIR"
          AS_VAR_COPY([arg_htmldir],[${arg_htmldir_name}])
          AS_IF([test "x${LALSUITE_BUILD}" = xtrue],[
            tagpath="\$(abs_top_builddir)/../${arg}/doxygen/out"
          ],[
            tagpath="${arg_htmldir}"
          ])
          DOXYGEN_TAGFILES="${DOXYGEN_TAGFILES} ${tagpath}/${arg}.tag=${tagpath}"
          DOXYGEN_INSTALL_DIRMAP="${DOXYGEN_INSTALL_DIRMAP} ${tagpath}:${arg_htmldir}"
        ]
      )
    done

    # configure MathJax
    AC_SUBST([DOXYGEN_MATHJAXDIR])
    AS_IF([test "x${DOXYGEN_MATHJAXDIR}" != x],[
      AS_IF([test ! -f "${DOXYGEN_MATHJAXDIR}/MathJax.js"],[
        AC_MSG_ERROR([no MathJax installation found in ${DOXYGEN_MATHJAXDIR}])
      ])
    ],[
      for dir in /usr/share/javascript/mathjax; do
        AC_MSG_CHECKING([for MathJax installation in ${dir}])
        AS_IF([test -f "${dir}/MathJax.js"],[
          AC_MSG_RESULT([yes])
          DOXYGEN_MATHJAXDIR="${dir}"
          break
        ],[
          AC_MSG_RESULT([no])
        ])
      done
      AS_IF([test "x${DOXYGEN_MATHJAXDIR}" = x],[
        DOXYGEN_MATHJAXDIR='https://cdn.mathjax.org/mathjax/latest'
        AC_MSG_NOTICE([using MathJax CDN at ${DOXYGEN_MATHJAXDIR}])
      ])
    ])

  ])
  # end $0
])

AC_DEFUN([LALSUITE_CHECK_SIMD],[
  # $0: check for SIMD extensions

  # list of recognised SIMD instruction sets
  m4_define([simd_isets],[m4_normalize([
    [SSE],[SSE2],[SSE3],[SSSE3],[SSE4.1],[SSE4.2],
    [AVX],[AVX2]
  ])])

  # push compiler environment
  LALSUITE_PUSH_UVARS
  LALSUITE_CLEAR_UVARS
  AC_LANG_PUSH([C])

  # check compiler support for each SIMD instruction set
  simd_supported=
  m4_foreach([iset],simd_isets,[
    m4_pushdef([option],m4_translit(iset,[A-Z.],[a-z.]))
    m4_pushdef([symbol],m4_translit(iset,[A-Z.],[A-Z_]))

    # assume supported until test fails, in which break out of loop
    for iset_supported in yes; do

      # check -m[]option flag
      LALSUITE_CHECK_COMPILE_FLAGS([-m[]option],[],[iset_supported=no])
      AS_IF([test x"${iset_supported}" = xno],[break])

      # check that compiler defines __]symbol[__ preprocessor symbol
      AC_MSG_CHECKING([whether ]_AC_LANG[ compiler defines __]symbol[__ with -m[]option])
      CFLAGS="${uvar_orig_prefix[]CFLAGS} -m[]option"
      AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM([],[[
#if !defined(__]]symbol[[__)
#error Preprocessor macro not defined by compiler
#endif
]])
      ],[
        AC_MSG_RESULT([yes])
      ],[
        AC_MSG_RESULT([no])
        iset_supported=no
      ])
      AS_IF([test x"${iset_supported}" = xno],[break])

      # check that compiler compiles floating-point math with -m[]option
      AC_MSG_CHECKING([whether ]_AC_LANG[ compiler compiles floating-point math with -m[]option])
      CFLAGS="${uvar_orig_prefix[]CFLAGS} -m[]option"
      AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM([
AC_INCLUDES_DEFAULT
#include <math.h>
],[[
double volatile a = 1.2;
double volatile b = 3.4;
double volatile c = a * b;
double volatile d = round(c);
]])
      ],[
        AC_MSG_RESULT([yes])
      ],[
        AC_MSG_RESULT([no])
        iset_supported=no
      ])
      AS_IF([test x"${iset_supported}" = xno],[break])

    done

    # define Automake conditional HAVE_<SIMD>_COMPILER
    AM_CONDITIONAL([HAVE_]symbol[_COMPILER],[test x"${iset_supported}" = xyes])
    AM_COND_IF([HAVE_]symbol[_COMPILER],[

      # define config.h preprocessor symbol HAVE_<SIMD>_COMPILER
      AC_DEFINE([HAVE_]symbol[_COMPILER],[1],[Define to 1 if compiler supports ]iset[ SIMD extensions])

      # substitute Automake variable <SIMD>_CFLAGS
      AC_SUBST(symbol[_CFLAGS],["-DSIMD_INSTRSET=]symbol[ -m[]option"])

      # add to list of supported instruction sets
      simd_supported="${simd_supported} iset"

      LALSUITE_ADD_TESTS_ENV_CONFIG_VAR([HAVE_]symbol[_COMPILER],[true])
    ],[
      LALSUITE_ADD_TESTS_ENV_CONFIG_VAR([HAVE_]symbol[_COMPILER],[false])
    ])

    m4_popdef([option])
    m4_popdef([symbol])
  ])

  # string listing all the SIMD extensions supported by the compiler
  simd_supported=`echo ${simd_supported}`
  AC_DEFINE_UNQUOTED([HAVE_SIMD_COMPILER],["${simd_supported}"],
    [Define to a string listing all the SIMD extensions supported by the ]_AC_LANG[ compiler]
  )
  AS_IF([test x"${simd_supported}" = x],[
    AC_MSG_NOTICE([]_AC_LANG[ compiler does not support SIMD instruction sets])
  ],[
    AC_MSG_NOTICE([]_AC_LANG[ compiler supports SIMD instruction sets: ${simd_supported}])
  ])

  # pop compiler environment
  AC_LANG_POP([C])
  LALSUITE_POP_UVARS

  # end $0
])

AC_DEFUN([LALSUITE_ENABLE_CFITSIO],[
  # $0: enable/disable cfitsio library
  AC_ARG_ENABLE(
    [cfitsio],
    AS_HELP_STRING([--enable-cfitsio],[compile code that requires cfitsio library [default=yes]]),
    AS_CASE(["${enableval}"],
      [yes],[cfitsio=true],
      [no],[cfitsio=false],
      AC_MSG_ERROR([bad value ${enableval} for --enable-cfitsio])
    ),
    [cfitsio=true]
  )
  # end $0
])

AC_DEFUN([LALSUITE_USE_CFITSIO],[
  # $0: check for cfitsio library
  AS_IF([test "x${cfitsio}" = xtrue],[
    LALSUITE_PUSH_UVARS
    PKG_CHECK_MODULES([CFITSIO],[cfitsio],[true],[false])
    LALSUITE_ADD_FLAGS([C],[${CFITSIO_CFLAGS}],[${CFITSIO_LIBS}])
    AC_CHECK_LIB([cfitsio],[ffopen],[],[cfitsio=false])
    AC_CHECK_HEADER([fitsio.h],[],[cfitsio=false])
    AC_CHECK_FUNCS([fffree])
    AC_CHECK_DECLS([fffree],[],[],[AC_INCLUDES_DEFAULT
#include <fitsio.h>
])
    LALSUITE_POP_UVARS
  ])
  AS_IF([test "x${cfitsio}" = xtrue],[
    LALSUITE_ADD_FLAGS([C],[${CFITSIO_CFLAGS}],[${CFITSIO_LIBS}])
  ])
  LALSUITE_ENABLE_MODULE([CFITSIO])
  # end $0
])

AC_DEFUN([LALSUITE_CHECK_PAGER],[
  # $0: check for pager programs and required functions
  AC_ARG_WITH(
    [pager],
    AS_HELP_STRING([--with-pager=PAGER],[specify pager program [default: less/more]]),
    [
      AS_CASE([${with_pager}],
        [no],[pager="false"],
        [yes],[pager="less -F -R -S -X more"],
        [pager="${with_pager}"]
      )
    ],[
      pager="less -F -R -S -X more"
    ]
  )
  AC_SUBST([PAGER],[])
  AC_SUBST([PAGER_CPPFLAGS],[])
  AS_IF([test "x${pager}" != xfalse],[
    for pager_arg in ${pager}; do
      AS_CASE([${pager_arg}],
        [-*],[
          AS_IF([test "x${PAGER}" != x],[
            AC_MSG_CHECKING([if option ${pager_arg} works with ${PAGER}])
            AS_IF([echo | ${PAGER} ${pager_arg} >/dev/null 2>&1],[
              AC_MSG_RESULT([yes])
              PAGER="${PAGER} ${pager_arg}"
            ],[
              AC_MSG_RESULT([no])
            ])
          ])
        ],[
          AS_IF([test "x${PAGER}" != x],[
            break
          ],[
            AC_PATH_PROGS([PAGER],[${pager_arg}],[])
          ])
        ]
      )
    done
  ])
  AS_IF([test "x${PAGER}" != x],[
      PAGER_CPPFLAGS="-DPAGER='\"\$(PAGER)\"'"
    AC_CHECK_FUNCS([popen pclose])
  ])
  # end $0
])

AC_DEFUN([LALSUITE_ENABLE_HELP2MAN],[
  # $0: check for help2man utility
  AC_PATH_PROG([HELP2MAN], [help2man])
  # Require a minimum version of help2man to support the
  # --version-string and --no-discard-stderr flags.
  AS_IF(
    [test -n "${HELP2MAN}"],
    [
      HELP2MAN_MIN_VERSION="1.37"
      AC_MSG_CHECKING([for help2man >= ${HELP2MAN_MIN_VERSION}])
      HELP2MAN_VERSION=`${HELP2MAN} --version --version | head -n 1 | cut -d ' ' -f 3`
      AX_COMPARE_VERSION(
        [${HELP2MAN_VERSION}], [ge], [${HELP2MAN_MIN_VERSION}],
        [AC_MSG_RESULT([yes])],
        [
          AC_MSG_RESULT([no])
          HELP2MAN=
        ]
      )
    ]
  )
  AC_ARG_ENABLE(
    [help2man],
    AS_HELP_STRING([--enable-help2man],[automatically generate man pages with help2man @<:@default=yes@:>@]),
    AS_CASE([${enableval}],
      [yes], [],
      [no], [HELP2MAN=],
      AC_MSG_ERROR([bad value ${enableval} for --enable-help2man])
    )
  )
  AC_SUBST([HELP2MAN], ["${HELP2MAN}"])
  AS_IF([test -n "${HELP2MAN}"], [help2man=true], [help2man=false])
  LALSUITE_ENABLE_MODULE([HELP2MAN])
  # end $0
])

AC_DEFUN([LALSUITE_ENABLE_OPENMP],[
  # $0: check for OpenMP support
  # check for OpenMP
  AC_OPENMP

  # check that compiler supports C99 variable length arrays in OpenMP for loops.
  # Apple's llvm-gcc-4.2 is buggy and does not.
  AS_IF([test "x$OPENMP_CFLAGS" != "x"],
    [
      AC_MSG_CHECKING([if compiler supports C99 VLAs in OpenMP for loops])
      LALSUITE_PUSH_UVARS
      CFLAGS="$OPENMP_CFLAGS"
      AC_LINK_IFELSE([
          AC_LANG_PROGRAM(, [
            int i, n = 10;
            #pragma omp parallel for
            for (i = 0; i < 10; i ++)
            { int a@<:@n@:>@; }
          ])
        ],
        [AC_MSG_RESULT([yes])],
        [
          AC_MSG_RESULT([no, disabling OpenMP])
          OPENMP_CFLAGS=
        ]
      )
      LALSUITE_POP_UVARS
    ]
  )

  # Disable OpenMP by default.
  # FIXME: OpenMP should be on by default, but it breaks condor_compiling lalapps.
  AS_IF(
    [test "x$enable_openmp" = "x"],
    [OPENMP_CFLAGS=]
  )
  # add compiler flags for OpenMP
  AS_IF(
    [test "x$OPENMP_CFLAGS" = "x" -a "x$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" != "xnone needed"],
    [
      openmp="false"
      LALSUITE_ADD_FLAGS([C],[-Wno-unknown-pragmas],[])
    ],
    [
      openmp="true"
      LALSUITE_ADD_FLAGS([C],[${OPENMP_CFLAGS}],[])
    ]
  )

  LALSUITE_ENABLE_MODULE([OPENMP])
  # end $0
])

AC_DEFUN([LALSUITE_ENABLE_MPI],
[AC_ARG_ENABLE(
  [mpi],
  AS_HELP_STRING([--enable-mpi],[compile using MPI for supported codes [default=no]]),
  [ case "${enableval}" in
      yes) mpi=true;;
      no)  mpi=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-mpi) ;;
    esac
  ], [ mpi=false ] )
])
