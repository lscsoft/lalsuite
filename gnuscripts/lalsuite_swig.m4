# -*- mode: autoconf; -*-
# lalsuite_swig.m4 - SWIG configuration
# Author: Karl Wette, 2011--2014
#
# serial 73

AC_DEFUN([_LALSUITE_CHECK_SWIG_VERSION],[
  # $0: check the version of $1, and store it in ${swig_version}
  swig_version=
  swig_version_output=[`$1 -version 2>/dev/null`]
  AS_IF([test $? -eq 0],[
    swig_version_regex=['s|^ *SWIG [Vv]ersion \([0-9.][0-9.]*\)|\1|p;d']
    swig_version=[`echo "${swig_version_output}" | ${SED} -e "${swig_version_regex}"`]
  ])
  AS_IF([test "x${swig_version}" = x],[
    AC_MSG_ERROR([could not determine version of $1])
  ])
  # end $0
])

AC_DEFUN([_LALSUITE_SWIG_CHECK_COMPILER_FLAGS],[
  # $0: check flags used to compile SWIG bindings
  LALSUITE_PUSH_UVARS
  LALSUITE_CLEAR_UVARS
  for flag in m4_normalize($2); do
    AC_MSG_CHECKING([if ]_AC_LANG[ compiler supports ${flag}])
    CFLAGS="-Werror ${flag}"
    CXXFLAGS="${CFLAGS}"
    AC_COMPILE_IFELSE([
      AC_LANG_PROGRAM([AC_INCLUDES_DEFAULT],[])
    ],[
      AC_MSG_RESULT([yes])
      $1="${$1} ${flag}"
    ],[
      AC_MSG_RESULT([no])
    ])
  done
  LALSUITE_POP_UVARS
  # end $0
])

AC_DEFUN([LALSUITE_ENABLE_SWIG],[
  # $0: enable SWIG bindings
  AC_REQUIRE([LALSUITE_CHECK_GIT_REPO])
  AM_COND_IF([HAVE_GIT_REPO],[swig_generate=true],[swig_generate=false])
  AC_ARG_ENABLE(
    [swig],
    AC_HELP_STRING(
      [--enable-swig],
      [generate SWIG bindings for all languages]
    ),[
      AS_CASE(["${enableval}"],
        [generate],[swig_build_all=true; swig_generate=true],
        [yes],[swig_build_all=true],
        [no],[swig_build_all=false],
        [AC_MSG_ERROR([invalid value "${enableval}" for --enable-swig])]
      )
    ],[
      swig_build_all=
    ]
  )
  swig_build_any=false
  LALSUITE_ENABLE_SWIG_LANGUAGE([Octave],[false],[LALSUITE_REQUIRE_CXX])
  LALSUITE_ENABLE_SWIG_LANGUAGE([Python],[false],[LALSUITE_REQUIRE_PYTHON([2.6])])
  AS_IF([test "${swig_generate}" = true],[
    # Python is required for running generate_swig_iface.py
    LALSUITE_REQUIRE_PYTHON([2.6])
    SWIG_GENERATE_ENABLE_VAL=ENABLED
  ],[
    SWIG_GENERATE_ENABLE_VAL=DISABLED
  ])
  # end $0
])

AC_DEFUN([LALSUITE_ENABLE_SWIG_LANGUAGE],[
  # $0: enable SWIG binding languages
  m4_pushdef([uppercase],m4_translit([$1],[a-z],[A-Z]))
  m4_pushdef([lowercase],m4_translit([$1],[A-Z],[a-z]))
  AC_ARG_ENABLE(
    [swig-]lowercase,
    AC_HELP_STRING(
      [--enable-swig-]lowercase,
      [generate SWIG binding for $1]
    ),[
      AS_CASE(["${enableval}"],
        [generate],[swig_build_]lowercase[=true; swig_generate=true],
        [yes],[swig_build_]lowercase[=true],
        [no],[swig_build_]lowercase[=false],
        [AC_MSG_ERROR([invalid value "${enableval}" for --enable-swig-]]lowercase[)]
      )
    ],[
      swig_build_]lowercase[=${swig_build_all:-$2}
    ]
  )
  AS_IF([test "${swig_build_]lowercase[}" = true],[
    swig_build_any=true
    SWIG_BUILD_]uppercase[_ENABLE_VAL=ENABLED
    $3
  ],[
    SWIG_BUILD_]uppercase[_ENABLE_VAL=DISABLED
  ])
  AM_CONDITIONAL([SWIG_BUILD_]uppercase,[test "${swig_build_]lowercase[}" = true])
  m4_popdef([uppercase])
  m4_popdef([lowercase])
  # end $0
])

AC_DEFUN([LALSUITE_USE_SWIG],[
  # $0: configure enabled SWIG bindings
  AS_IF([test "${swig_build_any}" = true],[

    # configure SWIG binding languages
    swig_min_version=2.0.11
    swig_min_version_info=""
    swig_src_files=""
    LALSUITE_USE_SWIG_OCTAVE
    LALSUITE_USE_SWIG_PYTHON

    # if SWIG bindings are being generated, check for SWIG binary
    # with version ${swig_min_version} or later; use value of
    # ${SWIG} if set, otherwise check common SWIG binary names
    AC_SUBST([SWIG])
    AS_IF([test "${swig_generate}" = true],[
      AS_IF([test "x${SWIG}" != x],[
        AC_MSG_CHECKING([if ${SWIG} version is at least ${swig_min_version}])
        _LALSUITE_CHECK_SWIG_VERSION([${SWIG}])
        LALSUITE_VERSION_COMPARE([${swig_version}],[<],[${swig_min_version}],[
          AC_MSG_RESULT([no (${swig_version})])
          AC_MSG_ERROR([SWIG version ${swig_min_version} or later is required ${swig_min_version_info}])
        ])
        AC_MSG_RESULT([yes (${swig_version})])
      ],[
        AC_PATH_PROGS_FEATURE_CHECK([SWIG],[swig swig2.0],[
          AC_MSG_CHECKING([if ${ac_path_SWIG} version is at least ${swig_min_version}])
          _LALSUITE_CHECK_SWIG_VERSION([${ac_path_SWIG}])
          LALSUITE_VERSION_COMPARE([${swig_version}],[>=],[${swig_min_version}],[
            ac_path_SWIG_found=true
            AC_MSG_RESULT([yes (${swig_version})])
            ac_cv_path_SWIG="${ac_path_SWIG}"
          ],[
            ac_path_SWIG_found=false
            AC_MSG_RESULT([no (${swig_version})])
          ])
        ],[
          AC_MSG_ERROR([SWIG version ${swig_min_version} or later is required ${swig_min_version_info}])
        ])
        SWIG="${ac_cv_path_SWIG}"
      ])
    ],[
      SWIG=false
    ])

    # if SWIG bindings are not being generated, check that the SWIG version
    # used to generate the bindings satisfies ${swig_min_version}
    src_file_swig_version_regex='s/^#.*define *SWIGVERSION *0x\([0-9][0-9]\)\([0-9][0-9]\)\([0-9][0-9]\).*$/\1 \2 \3/p'
    AS_IF([test "${swig_generate}" != true],[
      for file in ${swig_src_files}; do
        src_file="${srcdir}/swig/${file}"
        AS_IF([test -f "${src_file}"],[
          AC_MSG_CHECKING([if SWIG version ${swig_min_version} or later generated ${src_file}])
          src_file_swig_verargs=[`${SED} -n -e "${src_file_swig_version_regex}" "${src_file}"`]
          src_file_swig_version=[`printf '%d.%d.%d' ${src_file_swig_verargs}`]
          LALSUITE_VERSION_COMPARE([${src_file_swig_version}],[<],[${swig_min_version}],[
            AC_MSG_RESULT([no (${src_file_swig_version})])
            AC_MSG_ERROR([SWIG version ${swig_min_version} or later is required ${swig_min_version_info}])
          ])
          AC_MSG_RESULT([yes (${src_file_swig_version})])
        ],[
          AC_MSG_ERROR([${src_file} does not exist])
        ])
      done
    ])

    # extract -I and -D flags from LALSuite library preprocessor flags
    AC_SUBST([SWIG_CPPFLAGS],[])
    for flag in ${CPPFLAGS}; do
      AS_CASE([${flag}],
        [-I*|-D*],[SWIG_CPPFLAGS="${SWIG_CPPFLAGS} ${flag}"]
      )
    done

    # extract -L flags from LALSuite library linker flags
    AC_SUBST([SWIG_LDFLAGS],[])
    for flag in ${LDFLAGS}; do
      AS_CASE([${flag}],
        [-L*],[SWIG_LDFLAGS="${SWIG_LDFLAGS} ${flag}"]
      )
    done

    # determine various libtool parameters
    AC_SUBST([SWIG_LTLIBDIR],["${objdir}"])
    AC_SUBST([SWIG_SOEXT],[`module=yes; eval shrext=\"${shrext_cmds}\"; echo ${shrext}`])

    # substitute list of LALSuite bindings that these bindings depend on
    AC_MSG_CHECKING([for LALSuite binding dependencies])
    AC_SUBST([SWIG_DEPENDENCIES],[])
    for arg in ${lalsuite_libs}; do
      AS_CASE([${arg}],
        [lalsupport],[:],
        [SWIG_DEPENDENCIES="${SWIG_DEPENDENCIES} ${arg}"]
      )
    done
    AS_IF([test "x${SWIG_DEPENDENCIES}" = x],[
      AC_MSG_RESULT([none])
    ],[
      AC_MSG_RESULT([${SWIG_DEPENDENCIES}])
    ])

  ])
  AM_CONDITIONAL([SWIG_BUILD],[test "${swig_build_any}" = true])
  AM_CONDITIONAL([SWIG_GENERATE],[test "${swig_generate}" = true])
  # end $0
])

AC_DEFUN([LALSUITE_USE_SWIG_LANGUAGE],[
  # $0: configure SWIG binding languages
  m4_pushdef([uppercase],m4_translit([$1],[a-z],[A-Z]))
  m4_pushdef([lowercase],m4_translit([$1],[A-Z],[a-z]))
  AS_IF([test "${swig_build_]lowercase[}" = true],[
    swig_build=true
    swig_src_files="${swig_src_files} swiglal_[]lowercase[].$2"
    $3
  ])
  m4_popdef([uppercase])
  m4_popdef([lowercase])
  # end $0
])

AC_DEFUN([LALSUITE_USE_SWIG_OCTAVE],[
  # $0: configure SWIG Octave bindings
  LALSUITE_USE_SWIG_LANGUAGE([Octave],[cpp],[

    # check for Octave
    AC_PATH_PROG(OCTAVE,[octave],[],[])
    AS_IF([test "x${OCTAVE}" = x],[
      AC_MSG_ERROR([could not find octave in PATH])
    ])

    # check for Octave utilities octave-config and mkoctfile
    octave_dir=`AS_DIRNAME(["${OCTAVE}"])`
    AC_MSG_CHECKING([for octave-config])
    octave_cfg="${octave_dir}/octave-config"
    AS_IF([test -x "${octave_cfg}"],[
      AC_MSG_RESULT([${octave_cfg}])
    ],[
      AC_MSG_RESULT([not found])
      AC_MSG_ERROR([could not find octave-config in ${octave_dir}])
    ])
    octave_cfg="env - ${octave_cfg}"
    AC_MSG_CHECKING([for mkoctfile])
    mkoctfile="${octave_dir}/mkoctfile"
    AS_IF([test -x "${mkoctfile}"],[
      AC_MSG_RESULT([${mkoctfile}])
    ],[
      AC_MSG_RESULT([not found])
      AC_MSG_ERROR([could not find mkoctfile in ${octave_dir}])
    ])
    mkoctfile="env - ${mkoctfile}"

    # check Octave version
    octave_min_version=3.2.0
    AC_MSG_CHECKING([${OCTAVE} version])
    octave_version=[`${octave_cfg} -p VERSION 2>/dev/null`]
    AS_IF([test "x${octave_version}" = x],[
      AC_MSG_ERROR([could not determine ${OCTAVE} version])
    ])
    AC_MSG_RESULT([${octave_version}])
    LALSUITE_VERSION_COMPARE([${octave_version}],[<],[${octave_min_version}],[
      AC_MSG_ERROR([Octave version ${octave_min_version} or later is required])
    ])
    LALSUITE_VERSION_COMPARE([${octave_version}],[>=],[3.8.0],[
      swig_min_version=2.0.12
      swig_min_version_info="for Octave version ${octave_version}"
    ])

    # determine where to install Octave bindings: take versioned site .oct file directory given by octave-config,
    # and strip off prefix; thus, if LALSuite is installed in the same directory as Octave, .oct files will be
    # found by Octave without having to add to OCTAVE_PATH
    AC_MSG_CHECKING([${OCTAVE} .oct installation directory])
    octave_prefix=[`${octave_cfg} -p PREFIX 2>/dev/null | ${SED} -e 's|/*$||'`]
    octexecdir=[`${octave_cfg} -p LOCALVEROCTFILEDIR 2>/dev/null | ${SED} -e 's|/*$||'`]
    octexecdir=[`echo ${octexecdir} | ${SED} -e "s|^${octave_prefix}/||"`]
    AS_IF([test "x`echo ${octexecdir} | ${SED} -n -e '\|^/|p'`" != x],[
      AC_MSG_ERROR([could not build relative path from "${octexecdir}"])
    ])
    octexecdir='${prefix}'/"${octexecdir}"
    AC_MSG_RESULT([${octexecdir}])
    AC_SUBST([octexecdir])

    # check that wrappings are being compiled with the same C++ compiler used to compile Octave itself
    AC_MSG_CHECKING([C++ compiler used for building ${OCTAVE}])
    octave_CXX=`${mkoctfile} -p CXX 2>/dev/null`
    AC_MSG_RESULT([${octave_CXX}])
    CXX_version_regex=['s|([^)]*) *||g;s|^ *||g;s| *$||g']
    octave_CXX_version=`${octave_CXX} --version 2>/dev/null | ${SED} -n -e '1p' | ${SED} -e "${CXX_version_regex}"`
    _AS_ECHO_LOG([octave_CXX_version: '${octave_CXX_version}'])
    lalsuite_CXX_version=`${CXX} --version 2>/dev/null | ${SED} -n -e '1p' | ${SED} -e "${CXX_version_regex}"`
    _AS_ECHO_LOG([lalsuite_CXX_version: '${lalsuite_CXX_version}'])
    AS_IF([test "x${lalsuite_CXX_version}" != "x${octave_CXX_version}"],[
      AC_MSG_ERROR([configured C++ compiler "${CXX}" differs from ${OCTAVE} C++ compiler "${octave_CXX}"])
    ])

    # determine Octave preprocessor flags
    AC_SUBST([SWIG_OCTAVE_CPPFLAGS],[])
    AC_SUBST([SWIG_OCTAVE_CPPFLAGS_IOCTAVE],[])
    for arg in CPPFLAGS INCFLAGS; do
      for flag in `${mkoctfile} -p ${arg} 2>/dev/null`; do
        AS_CASE([${flag}],
          [-I*/octave],[SWIG_OCTAVE_CPPFLAGS_IOCTAVE="${flag}"],
          [SWIG_OCTAVE_CPPFLAGS="${SWIG_OCTAVE_CPPFLAGS} ${flag}"]
        )
      done
    done

    # determine Octave compiler flags
    AC_SUBST([SWIG_OCTAVE_CXXFLAGS],[])
    for arg in CXXPICFLAG ALL_CXXFLAGS; do
      SWIG_OCTAVE_CXXFLAGS="${SWIG_OCTAVE_CXXFLAGS} "`${mkoctfile} -p ${arg} 2>/dev/null`
    done
    AC_LANG_PUSH([C++])
    _LALSUITE_SWIG_CHECK_COMPILER_FLAGS([SWIG_OCTAVE_CXXFLAGS],[
      -Wno-uninitialized -Wno-unused-variable -Wno-unused-but-set-variable
      -Wno-tautological-compare -fno-strict-aliasing -g
      -O0 -Wp[,]-U_FORTIFY_SOURCE
    ])
    AC_LANG_POP([C++])

    # determine Octave linker flags
    AC_SUBST([SWIG_OCTAVE_LDFLAGS],[])
    for arg in LFLAGS LIBOCTINTERP LIBOCTAVE LIBCRUFT OCT_LINK_OPTS OCT_LINK_DEPS; do
      SWIG_OCTAVE_LDFLAGS="${SWIG_OCTAVE_LDFLAGS} "`${mkoctfile} -p ${arg} 2>/dev/null`
    done

    # check for Octave headers
    AC_LANG_PUSH([C++])
    LALSUITE_PUSH_UVARS
    CPPFLAGS="${SWIG_OCTAVE_CPPFLAGS_IOCTAVE} ${SWIG_OCTAVE_CPPFLAGS}"
    AC_CHECK_HEADERS([octave/oct.h],[],[
      AC_MSG_ERROR([could not find the header "octave/oct.h"])
    ],[
      AC_INCLUDES_DEFAULT
    ])
    LALSUITE_POP_UVARS
    AC_LANG_POP([C++])

  ])
  # end $0
])

AC_DEFUN([LALSUITE_USE_SWIG_PYTHON],[
  # $0: configure SWIG Python bindings
  LALSUITE_USE_SWIG_LANGUAGE([Python],[c],[

    # check for distutils
    AC_MSG_CHECKING([for distutils])
    cat <<EOD | ${PYTHON} - 2>/dev/null
import distutils
EOD
    AS_IF([test $? -ne 0],[
      AC_MSG_ERROR([could not import distutils])
    ])
    AC_MSG_RESULT([yes])

    # check for NumPy
    numpy_min_version=1.3
    AC_MSG_CHECKING([for NumPy])
    numpy_version=[`cat <<EOD | ${PYTHON} - 2>/dev/null
import numpy
print(numpy.__version__)
EOD`]
    AS_IF([test $? -ne 0],[
      AC_MSG_ERROR([could not import NumPy])
    ])
    AC_MSG_RESULT([yes])

    # check NumPy version
    AC_MSG_CHECKING([NumPy version])
    LALSUITE_VERSION_COMPARE([${numpy_version}],[<],[${numpy_min_version}],[
      AC_MSG_ERROR([NumPy version ${numpy_min_version} or later is required])
    ])
    AC_MSG_RESULT([${numpy_version}])

    # determine Python preprocessor flags
    AC_SUBST([SWIG_PYTHON_CPPFLAGS],[])
    python_out=[`cat <<EOD | ${PYTHON} - 2>/dev/null
import sys
import distutils.sysconfig as cfg
import numpy.lib.utils as npyutil
sys.stdout.write( '-I' + cfg.get_python_inc())
sys.stdout.write(' -I' + cfg.get_python_inc(plat_specific=1))
sys.stdout.write(' -I' + npyutil.get_include())
EOD`]
    AS_IF([test $? -ne 0],[
      AC_MSG_ERROR([could not determine Python preprocessor flags])
    ])
    SWIG_PYTHON_CPPFLAGS="${SWIG_PYTHON_CPPFLAGS} ${python_out}"

    # determine Python compiler flags
    AC_SUBST([SWIG_PYTHON_CFLAGS],[])
    python_out=[`cat <<EOD | ${PYTHON} - 2>/dev/null
import sys
import distutils.sysconfig as cfg
cflags = cfg.get_config_var('CFLAGS').split()
cflags = [f for f in cflags if f != '-DNDEBUG']
sys.stdout.write(" ".join(cflags))
EOD`]
    AS_IF([test $? -ne 0],[
      AC_MSG_ERROR([could not determine Python compiler flags])
    ])
    SWIG_PYTHON_CFLAGS="${SWIG_PYTHON_CFLAGS} ${python_out}"
    AC_LANG_PUSH([C])
    _LALSUITE_SWIG_CHECK_COMPILER_FLAGS([SWIG_PYTHON_CFLAGS],[
      -Wno-uninitialized -Wno-unused-variable -Wno-unused-but-set-variable
      -Wno-tautological-compare -fno-strict-aliasing -g
    ])
    AC_LANG_POP([C])

    # determine Python linker flags
    AC_SUBST([SWIG_PYTHON_LDFLAGS],[])
    python_out=[`cat <<EOD | ${PYTHON} - 2>/dev/null
import sys, os
import distutils.sysconfig as cfg
sys.stdout.write(cfg.get_config_var('LINKFORSHARED'))
sys.stdout.write(' -L' + cfg.get_python_lib())
sys.stdout.write(' -L' + cfg.get_python_lib(plat_specific=1))
sys.stdout.write(' -L' + cfg.get_python_lib(plat_specific=1,standard_lib=1))
sys.stdout.write(' -L' + cfg.get_config_var('LIBDIR'))
EOD`]
    AS_IF([test $? -ne 0],[
      AC_MSG_ERROR([could not determine Python linker flags])
    ])
    SWIG_PYTHON_LDFLAGS="${SWIG_PYTHON_LDFLAGS} ${python_out}"

    # check for Python and NumPy headers
    AC_LANG_PUSH([C])
    LALSUITE_PUSH_UVARS
    CPPFLAGS="${SWIG_PYTHON_CPPFLAGS}"
    AC_CHECK_HEADERS([Python.h],[],[
      AC_MSG_ERROR([could not find the header "Python.h"])
    ],[
      AC_INCLUDES_DEFAULT
    ])
    AC_CHECK_HEADERS([numpy/arrayobject.h],[],[
      AC_MSG_ERROR([could not find the header "numpy/arrayobject.h"])
    ],[
      AC_INCLUDES_DEFAULT
      #include <Python.h>
    ])
    LALSUITE_POP_UVARS
    AC_LANG_POP([C])

    # remove deprecated code in NumPy API >= 1.7
    SWIG_PYTHON_CPPFLAGS="${SWIG_PYTHON_CPPFLAGS} -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION"

    # check for declarations which may need compatibility code for NumPy API < 1.7
    AC_LANG_PUSH([C])
    LALSUITE_PUSH_UVARS
    CPPFLAGS="${SWIG_PYTHON_CPPFLAGS}"
    AC_CHECK_DECLS([NPY_ARRAY_WRITEABLE,PyArray_SetBaseObject],,,[
      AC_INCLUDES_DEFAULT
      #include <Python.h>
      #include <numpy/arrayobject.h>
    ])
    LALSUITE_POP_UVARS
    AC_LANG_POP([C])

  ])
  # end $0
])
