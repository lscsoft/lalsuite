# SWIG configuration
# Author: Karl Wette, 2011
#
# serial 14

# basic version string comparison
# can only handle numeric versions separated by periods
AC_DEFUN([LALSUITE_VERSION_COMPARE],[
  vcmp_awkprog='{n = 0; while(++n <= NF) { if($n > 99) printf "ERROR"; printf "%02u", $n; } }'
  vcmp_v1=[`echo $1 | ${SED} '/[^0-9.]/d' | ${AWK} -F . "${vcmp_awkprog}" | ${SED} '/ERROR/d'`]
  vcmp_v2=[`echo $2 | ${SED} '/[^0-9.]/d' | ${AWK} -F . "${vcmp_awkprog}" | ${SED} '/ERROR/d'`]
  AS_IF([test x${vcmp_v1} = x],[AC_MSG_ERROR([could not parse version string '$1'])])
  AS_IF([test x${vcmp_v2} = x],[AC_MSG_ERROR([could not parse version string '$2'])])
  AS_IF([test ${vcmp_v1} -lt ${vcmp_v2}],[$3])
  AS_IF([test ${vcmp_v1} -eq ${vcmp_v2}],[$4])
  AS_IF([test ${vcmp_v1} -gt ${vcmp_v2}],[$5])
])

# workaround to check whether SWIG modules are going to be
# built (and therefore a C++ compiler is required) before
# calling LALSUITE_PROG_CC_CXX, since LALSUITE_ENABLE_SWIG
# must appear after LALSUITE_PROG_CC_CXX. to be fixed!
AC_DEFUN([LALSUITE_SWIG_REQUIRE_CXX],[
  AS_IF([test "${enable_swig}" = yes],[LALSUITE_REQUIRE_CXX])
  AS_IF([test "${enable_swig_octave}" = yes],[LALSUITE_REQUIRE_CXX])
  AS_IF([test "${enable_swig_python}" = yes],[LALSUITE_REQUIRE_CXX])
])

# SWIG setup and configuration
AC_DEFUN([LALSUITE_ENABLE_SWIG],[

  # minimum required SWIG version
  SWIG_MIN_VERSION=1.3.40

  # save and clear CPPFLAGS and LIBS
  swig_CPPFLAGS=${CPPFLAGS}
  swig_LIBS=${LIBS}
  CPPFLAGS=
  LIBS=

  # check for sed and awk
  AC_PROG_SED
  AC_PROG_AWK

  # check for MKDIR_P
  m4_ifdef([AC_PROG_MKDIR_P],[],[
    MKDIR_P='$(INSTALL) -d'
    AC_SUBST(MKDIR_P)
  ])

  # command line option to enable/disable all languages
  AC_ARG_ENABLE(
    [swig],
    AC_HELP_STRING(
      [--enable-swig],
      [generate SWIG wrappings for all languages]
    ),[
      case "${enableval}" in
        yes) swig_build_all=true;;
        no)  swig_build_all=false;;
        *)   AC_MSG_ERROR([invalid value '${enableval}' for --enable-swig]);;
      esac
    ],[
      swig_build_all=
    ]
  )

  # command line option to use specific SWIG binary
  AC_ARG_WITH(
    [swig],
    AC_HELP_STRING(
      [--with-swig],
      [specify SWIG binary (default: search \$PATH)]
    ),[
      AS_IF([test -f "${withval}"],[
        SWIG="${withval}"
      ],[
        AC_MSG_ERROR([file '${withval}' not found])
      ])
    ],[
      SWIG=
    ]
  )

  # are we are binding LAL itself, or one of the other LAL libraries?
  AS_IF([test x${PACKAGE_NAME} = xlal],[
    swig_is_lal=true
  ],[
    swig_is_lal=false
  ])

  # common SWIG interface headers (with LAL only)
  AS_IF([test ${swig_is_lal} = true],[
    SWIG_HEADERS=
    AC_SUBST(SWIG_HEADERS)
  ])

  # string to add to user environment setup scripts
  SWIG_USER_ENV=

  # configure SWIG target scripting languages
  swig_build=false
  LALSUITE_SWIG_LANGUAGES

  # check if any language was configured
  AM_CONDITIONAL(SWIG_BUILD,[test ${swig_build} = true])
  AS_IF([test ${swig_build} = true],[

    # check for swig binary
    AC_MSG_CHECKING([for swig])
    AS_IF([test "x${SWIG}" = x],[
      AC_PATH_PROGS(SWIG,[swig],[])
      AS_IF([test "x${SWIG}" = x],[
        AC_MSG_ERROR([could not find 'swig' in path])
      ])
    ])
    AC_MSG_RESULT([${SWIG}])

    # check for swig version
    AC_MSG_CHECKING([for swig version])
    swig_regex=['s|^ *SWIG [Vv]ersion \([0-9.][0-9.]*\) *$|\1|p;d']
    swig_version=[`${SWIG} -version | ${SED} "${swig_regex}"`]
    AS_IF([test "x${swig_version}" = x],[
      AC_MSG_ERROR([could not determine swig version])
    ])
    AC_MSG_RESULT([${swig_version}])

    # check if swig version is newer than required
    LALSUITE_VERSION_COMPARE([${SWIG_MIN_VERSION}],[${swig_version}],[],[],[
      AC_MSG_ERROR([require swig version >= ${SWIG_MIN_VERSION}])
    ])

    # check for perl binary
    AC_PATH_PROGS(PERL,[perl],[])
    AS_IF([test "x${PERL}" = x],[
      AC_MSG_ERROR([could not find 'perl' in path])
    ])

    # symbols to define when generating SWIG wrapping code
    SWIG_SWIG_DEFINES=
    AC_SUBST(SWIG_SWIG_DEFINES)

    # symbols to define when compiling SWIG wrapping code
    SWIG_CXX_DEFINES=
    AC_SUBST(SWIG_CXX_DEFINES)

    # are we are binding LAL itself, or one of the other LAL libraries?
    AS_IF([test ${swig_is_lal} = true],[
      SWIG_SWIG_DEFINES="${SWIG_SWIG_DEFINES} SWIGLAL_IS_LAL"
    ])

    # are we (not) in debugging mode?
    AS_IF([test x${enable_debug} = xno],[
      SWIG_SWIG_DEFINES="${SWIG_SWIG_DEFINES} SWIGLAL_NDEBUG"
      SWIG_CXX_DEFINES="${SWIG_CXX_DEFINES} SWIGLAL_NDEBUG"
    ])

    # try to figure out the underlying type of int64_t
    AC_CHECK_HEADERS([stdint.h],[],[
      AC_MSG_ERROR([could not find 'stdint.h'])
    ])
    AC_MSG_CHECKING([underlying type of int64_t])
    AC_LANG_PUSH([C++])
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM([AC_INCLUDES_DEFAULT],[
          int64_t i64 = 0; const long int *pli = &i64;
        ])
      ],[
        AC_MSG_RESULT([long int])
        swig_wordsize=SWIGWORDSIZE64
      ],[
        AC_COMPILE_IFELSE(
          [
            AC_LANG_PROGRAM([AC_INCLUDES_DEFAULT],[
              int64_t i64 = 0; const long long int *plli = &i64;
            ])
          ],[
            AC_MSG_RESULT([long long int])
            swig_wordsize=
          ],[
            AC_MSG_ERROR([could not determine underlying type of int64_t])
          ]
        )
      ]
    )
    AC_LANG_POP([C++])
    SWIG_SWIG_DEFINES="${SWIG_SWIG_DEFINES} ${swig_wordsize}"

    # ensure that all LAL library modules share type information
    SWIG_SWIG_DEFINES="${SWIG_SWIG_DEFINES} SWIG_TYPE_TABLE=swiglaltypetable"

    # make SWIG use C++ casts
    SWIG_SWIG_DEFINES="${SWIG_SWIG_DEFINES} SWIG_CPLUSPLUS_CAST"

    # define C99 constant and limit macros
    SWIG_CXX_DEFINES="${SWIG_CXX_DEFINES} __STDC_CONSTANT_MACROS __STDC_LIMIT_MACROS"               

    # common SWIG interface headers (with LAL only)
    AS_IF([test ${swig_is_lal} = true],[
      SWIG_HEADERS="${SWIG_HEADERS} \$(swig_srcdir)/swiglal-common.i"
      SWIG_HEADERS="${SWIG_HEADERS} \$(swig_srcdir)/swiglal-gsl.i"
      SWIG_HEADERS="${SWIG_HEADERS} \$(swig_srcdir)/swiglal-test.i"
    ])

    # string to add to user environment setup scripts
    AC_SUBST(SWIG_USER_ENV)

    # path SWIG should look in for header files:
    #  - keep any -I options in CPPFLAGS, without the -I prefix
    SWIG_INCLPATH=[`for n in ${swig_CPPFLAGS}; do echo $n | ${SED} 's|^-I||p;d'; done`]
    SWIG_INCLPATH=[`echo ${SWIG_INCLPATH}`]   # get rid of newlines
    AC_SUBST(SWIG_INCLPATH)

    # path SWIG should look in for (pre-installed) libraries:
    #  - keep any -L options in _LIB variables, without the -L prefix
    #  - keep any "lib*.la" files, replace filename with $objdir (pre-install)
    swig_all_libs="${LAL_LIBS} ${LALSUPPORT_LIBS} ${swig_LIBS}"
    SWIG_LIBPATH=[`for n in ${swig_all_libs}; do echo $n | ${SED} 's|^-L||p;d'; done`]
    SWIG_LIBPATH=[`echo ${SWIG_LIBPATH}`]   # get rid of newlines
    AC_SUBST(SWIG_LIBPATH)
    SWIG_PREINST_LIBPATH=[`for n in ${swig_all_libs}; do echo $n | ${SED} 's|lib[^/][^/]*\.la|'"${objdir}"'|p;d'; done`]
    SWIG_PREINST_LIBPATH=[`echo ${SWIG_PREINST_LIBPATH}`]   # get rid of newlines
    SWIG_PREINST_LIBPATH="${SWIG_PREINST_LIBPATH} \$(top_builddir)/lib/${objdir}"
    SWIG_PREINST_LIBPATH="${SWIG_PREINST_LIBPATH} \$(top_builddir)/src/${objdir}"
    SWIG_PREINST_LIBPATH="${SWIG_PREINST_LIBPATH} \$(top_builddir)/packages/support/src/${objdir}"
    AC_SUBST(SWIG_PREINST_LIBPATH)

    # deduce library load path to use when running check scripts prior to installation
    AS_IF([test ${build_vendor} = apple],[
      SWIG_LD_LIBPATH_NAME=DYLD_FALLBACK_LIBRARY_PATH
    ],[
      SWIG_LD_LIBPATH_NAME=LD_LIBRARY_PATH
    ])
    AC_SUBST(SWIG_LD_LIBPATH_NAME)

  ],[

    # if no SWIG languages were found
    SWIG_WRAPPINGS="NONE"

  ])

  # restore CPPFLAGS and LIBS
  CPPFLAGS=${swig_CPPFLAGS}
  LIBS=${swig_LIBS}

])

# tell the SWIG wrappings to use some feature
AC_DEFUN([LALSUITE_SWIG_USE],[
  SWIG_SWIG_DEFINES="${SWIG_SWIG_DEFINES} SWIGLAL_USE_$1"
])

# SWIG language configuration
AC_DEFUN([LALSUITE_SWIG_LANGUAGE],[

  # uppercase and lowercase language name
  m4_pushdef([uppercase],translit([$1],[a-z],[A-Z]))
  m4_pushdef([lowercase],translit([$1],[A-Z],[a-z]))

  # command line option to enable/disable $1
  AC_ARG_ENABLE(
    [swig-]lowercase,
    AC_HELP_STRING(
      [--enable-swig-]lowercase,
      [generate SWIG wrappings for $1]
    ),[
      case "${enableval}" in
        yes) swig_build_]lowercase[=true;;
        no)  swig_build_]lowercase[=false;;
        *)   AC_MSG_ERROR([invalid value '${enableval}' for --enable-swig-]]lowercase[);;
      esac
    ],[
      swig_build_]lowercase[=${swig_build_all:-false}
    ]
  )

  # check whether to configure $1
  AM_CONDITIONAL(SWIG_BUILD_[]uppercase,[test ${swig_build_]lowercase[} = true])
  AS_IF([test ${swig_build_]lowercase[} = true],[
    
    # at least one language was configured
    swig_build=true

    # set message string to indicate language will be built
    SWIG_]uppercase[_ENABLE_VAL=ENABLED

    # language-specific SWIG interface headers (with LAL only)
    AS_IF([test ${swig_is_lal} = true],[
      SWIG_]uppercase[_HEADERS="\$(swig_srcdir)/]lowercase[/swiglal-]lowercase[.i"
      AC_SUBST(SWIG_]uppercase[_HEADERS)
    ])

    # configure $1
    $2
    # $1 has been configured

  ],[
    SWIG_]uppercase[_ENABLE_VAL=DISABLED
  ])

  # clear M4 definitions
  m4_popdef([uppercase])
  m4_popdef([lowercase])

])

# SWIG languages
AC_DEFUN([LALSUITE_SWIG_LANGUAGES],[
  LALSUITE_SWIG_LANGUAGE_OCTAVE
  LALSUITE_SWIG_LANGUAGE_PYTHON
])

# SWIG octave configuration
AC_DEFUN([LALSUITE_SWIG_LANGUAGE_OCTAVE],[
  LALSUITE_SWIG_LANGUAGE([Octave],[

    # minimum required octave version
    OCTAVE_MIN_VERSION=3.2.0

    # check for octave-config binary
    AC_PATH_PROGS(OCTAVE_CONFIG,[octave-config],[])
    AS_IF([test "x${OCTAVE_CONFIG}" = x],[
      AC_MSG_ERROR([could not find 'octave-config' in path])
    ])

    # check for corresponding octave binary
    AC_MSG_CHECKING([for octave])
    OCTAVE=`${OCTAVE_CONFIG} -p BINDIR`/octave
    AS_IF([test -f "${OCTAVE}" && test -x "${OCTAVE}"],[],[
      AC_MSG_ERROR([could not find 'octave' in path])
    ])
    AC_MSG_RESULT([${OCTAVE}])
    AC_SUBST(OCTAVE)
    # add flags for silence and environment-independence
    OCTAVE="${OCTAVE} -qfH"

    # check for octave version
    AC_MSG_CHECKING([for octave version])
    octave_version=`${OCTAVE_CONFIG} --version`
    AS_IF([test "x${octave_version}" = x],[
      AC_MSG_ERROR([could not determine octave version])
    ])
    AC_MSG_RESULT([${octave_version}])

    # check if octave version is newer than required
    LALSUITE_VERSION_COMPARE([${OCTAVE_MIN_VERSION}],[${octave_version}],[],[],[
      AC_MSG_ERROR([require octave version >= ${OCTAVE_MIN_VERSION}])
    ])

    # determine where to install .oct files:
    # take site .oct install dir given by octave-config,
    # and strip off prefix; thus, if LAL is installed in
    # the same directory as octave, .oct files will be
    # found by octave without having to add to OCTAVE_PATH
    octave_prefix=[`${OCTAVE_CONFIG} -p PREFIX | ${SED} 's|/*$||'`]
    AC_MSG_CHECKING([for octave .oct installation directory])
    octave_localoctfiledir=[`${OCTAVE_CONFIG} -p LOCALOCTFILEDIR | ${SED} 's|/*$||'`]
    octave_octfiledir=[`echo ${octave_localoctfiledir} | ${SED} "s|^${octave_prefix}/||"`]
    AS_IF([test -n "`echo ${octave_octfiledir} | ${SED} -n '\|^/|p'`"],[
      AC_MSG_ERROR([could not build relative path from '${octave_octfiledir}'])
    ])
    AC_MSG_RESULT([\${prefix}/${octave_octfiledir}])
    AC_SUBST(octfiledir, [${prefix}/${octave_octfiledir}])

    # string to add to user environment setup scripts
    SWIG_USER_ENV="${SWIG_USER_ENV}"'prepend OCTAVE_PATH $(octfiledir)~E~O~L~'

  ])
])

# SWIG python configuration
AC_DEFUN([LALSUITE_SWIG_LANGUAGE_PYTHON],[
  LALSUITE_SWIG_LANGUAGE([Python],[

    # check for python
    AM_PATH_PYTHON([2.4])

    # check for numpy
    AC_MSG_CHECKING([for numpy])
    ${PYTHON} -c "import numpy" 2>/dev/null
    AS_IF([test $? -ne 0],[
      AC_MSG_ERROR([could not import numpy])
    ])
    AC_MSG_RESULT([yes])

    # string to add to user environment setup scripts
    SWIG_USER_ENV="${SWIG_USER_ENV}"'prepend PYTHONPATH $(pyexecdir)~E~O~L~'

  ])
])
