AC_PREREQ([2.63])
AC_INIT([LALSuite],[7.25.1],[lal-discuss@ligo.org])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_PREFIX_DEFAULT(["`pwd`/_inst"])
AC_CONFIG_FILES([ \
  Makefile \
  doxygen/Makefile \
  wheel/Makefile \
  wheel/setup.py \
])
AM_INIT_AUTOMAKE([ \
  1.11 \
  foreign \
  subdir-objects \
  color-tests \
  parallel-tests \
  tar-pax \
  dist-xz \
  no-dist-gzip \
])

LALSUITE_DISTCHECK_CONFIGURE_FLAGS

# use silent build rules
AM_SILENT_RULES([yes])

# check for programs
AC_PROG_AWK
AC_PROG_SED

LALSUITE_ENABLE_NIGHTLY
LALSUITE_CHECK_PYTHON([3.5])
LALSUITE_CHECK_GIT_REPO

# provide LAL library enable/disable options
LALSUITE_ENABLE_ALL_LAL
lal=true
LALSUITE_ENABLE_LALFRAME
LALSUITE_ENABLE_LALMETAIO
LALSUITE_ENABLE_LALSIMULATION
LALSUITE_ENABLE_LALBURST
LALSUITE_ENABLE_LALINSPIRAL
LALSUITE_ENABLE_LALINFERENCE
LALSUITE_ENABLE_LALPULSAR
LALSUITE_ENABLE_LALAPPS

# enable Doxygen top-level documentation
LALSUITE_ENABLE_DOXYGEN

# directories of top-level build and source
lalsuite_abs_top_builddir="`pwd`"
lalsuite_abs_top_srcdir='$(abs_top_srcdir)/..'

# indicate to configure.ac in subdirectories that this is a top-level build
LALSUITE_BUILD="true"
export LALSUITE_BUILD

# list of LALSuite packages being built
LALSUITE_PACKAGES=
export LALSUITE_PACKAGES

# configure a lalsuite package
lalsuite_libs=
lalsuite_package_configure_deps=
AC_DEFUN([lalsuite_config_subdir],[
  m4_pushdef([lowercase],m4_translit([[$1]], [A-Z], [a-z]))
  m4_pushdef([uppercase],m4_translit([[$1]], [a-z], [A-Z]))

  # if $1 is enabled
  AS_IF([test "x${lowercase}" = xtrue],[

    # add to list of packages being built
    LALSUITE_PACKAGES="${LALSUITE_PACKAGES} $1"

    # check that environment for enabled dependent pacakges has been set
    m4_foreach_w([dep],$2,[
      m4_pushdef([dep_lowercase],m4_translit(dep, [A-Z], [a-z]))
      m4_pushdef([dep_uppercase],m4_translit(dep, [a-z], [A-Z]))
      AS_IF([test "x${dep_lowercase}" = xtrue && test "x${dep_uppercase[]_[]CFLAGS}" = x],[
        AC_MSG_ERROR([$1 depends on dep but dep environment is not set; [lalsuite_config_subdir()] calls may be in wrong dependency order])
      ])
      m4_popdef([dep_lowercase])
      m4_popdef([dep_uppercase])
    ])

    # re-run top-level ./00boot if any package-level ./configure.ac have changed
    # - this is needed to pick up e.g. any new inter-package dependencies
    lalsuite_package_configure_deps="${lalsuite_package_configure_deps} "'$(top_srcdir)/lowercase/configure.ac'

    # re-build AUTHORS if any package-level AUTHORS have changed
    lalsuite_package_AUTHORS_deps="${lalsuite_package_AUTHORS_deps} "'$(top_srcdir)/lowercase/AUTHORS'

    m4_if($1,[LALApps],[],[

      # add to list of libraries being built
      lalsuite_libs="${lalsuite_libs} lowercase"

      # export precious environment variables
      uppercase[]_LIBS="${lalsuite_abs_top_builddir}/lowercase/lib/lib[]lowercase.la"
      uppercase[]_CFLAGS="-I${lalsuite_abs_top_builddir}/lowercase/include"
      uppercase[]_HTMLDIR="${htmldir}/../lowercase"
      m4_foreach_w([dep],$2,[
        m4_pushdef([dep_uppercase],m4_translit(dep, [a-z], [A-Z]))
        for flag in ${dep_uppercase[]_CFLAGS}; do
          AS_CASE([${flag}],
            [-I*],[
              AS_CASE([" ${uppercase[]_CFLAGS} "],
                [*" ${flag} "*],[:],
                [uppercase[]_CFLAGS="${uppercase[]_CFLAGS} ${flag}"]
              )
            ]
          )
        done
        m4_popdef([dep_uppercase])
      ])
      export uppercase[]_LIBS uppercase[]_CFLAGS uppercase[]_HTMLDIR

      # export precious environment paths
      uppercase[]_BIN_PATH="${lalsuite_abs_top_builddir}/lowercase/bin"
      uppercase[]_DATA_PATH="${lalsuite_abs_top_srcdir}/lowercase/lib:${lalsuite_abs_top_srcdir}/lowercase/test"
      uppercase[]_OCTAVE_PATH="${lalsuite_abs_top_builddir}/lowercase/octave"
      uppercase[]_PYTHON_PATH="${lalsuite_abs_top_builddir}/lowercase/python"
      m4_foreach_w([suffix],DATA_PATH OCTAVE_PATH PYTHON_PATH,[
        m4_foreach_w([dep],$2,[
          m4_pushdef([dep_uppercase],m4_translit(dep, [a-z], [A-Z]))
          LALSUITE_ADD_PATH(uppercase[]_[]suffix,"${dep_uppercase[]_[]suffix}")
          m4_popdef([dep_uppercase])
        ])
      ])
      export uppercase[]_BIN_PATH uppercase[]_DATA_PATH uppercase[]_OCTAVE_PATH uppercase[]_PYTHON_PATH

    ])

    # configure $1
    AC_CONFIG_SUBDIRS(lowercase)

    # set enable string
    uppercase[]_ENABLE_VAL=ENABLED
  ],[
    # set disable string
    uppercase[]_ENABLE_VAL=DISABLED
  ])

  m4_popdef([lowercase])
  m4_popdef([uppercase])
])

# lal and lalsupport are always configured
lalsuite_config_subdir([LAL],[])
LALSUPPORT_LIBS="${lalsuite_abs_top_builddir}/lal/lib/support/liblalsupport.la"
LALSUPPORT_CFLAGS="-I${lalsuite_abs_top_builddir}/lal/include"
LALSUPPORT_BIN_PATH="${LAL_BIN_PATH}"
LALSUPPORT_DATA_PATH="${LAL_DATA_PATH}"
LALSUPPORT_OCTAVE_PATH="${LAL_OCTAVE_PATH}"
LALSUPPORT_PYTHON_PATH="${LAL_PYTHON_PATH}"
LALSUPPORT_HTMLDIR="${htmldir}/../lal"
export LALSUPPORT_LIBS LALSUPPORT_CFLAGS LALSUPPORT_HTMLDIR
export LALSUPPORT_BIN_PATH LALSUPPORT_DATA_PATH LALSUPPORT_OCTAVE_PATH LALSUPPORT_PYTHON_PATH

# configure optional packages
# - 2nd argument lists library dependencies
#   * used by lalsuite_config_subdir() to generate DATA/OCTAVE/PYTHON_PATHs including dependencies
#   * used by gnuscripts/check_library_dependencies to check for consistent dependencies throughout build system
lalsuite_config_subdir([LALFrame],[LAL])
lalsuite_config_subdir([LALMetaIO],[LAL])
lalsuite_config_subdir([LALSimulation],[LAL])
lalsuite_config_subdir([LALBurst],[LAL LALMetaIO LALSimulation])
lalsuite_config_subdir([LALInspiral],[LAL LALFrame LALMetaIO LALSimulation LALBurst])
lalsuite_config_subdir([LALInference],[LAL LALFrame LALMetaIO LALSimulation LALBurst LALInspiral])
lalsuite_config_subdir([LALPulsar],[LAL LALFrame LALSimulation LALInference])
lalsuite_config_subdir([LALApps],[LAL LALFrame LALMetaIO LALSimulation LALBurst LALInspiral LALInference LALPulsar])

# configure Doxygen top-level documentation
# - needs to come after all LALSuite libraries have been configured
LALSUITE_USE_DOXYGEN
DOXYGEN_TAGFILES=`echo "${DOXYGEN_TAGFILES}" | ${SED} 's|\$(abs_top_builddir)/\.\./|$(abs_top_builddir)/|g'`
DOXYGEN_INSTALL_DIRMAP=`echo "${DOXYGEN_INSTALL_DIRMAP}" | ${SED} 's|\$(abs_top_builddir)/\.\./|$(abs_top_builddir)/|g'`

AC_SUBST([LALSUITE_PACKAGES])
AC_SUBST([ac_configure_args])
AC_SUBST([lalsuite_package_configure_deps])
AC_SUBST([lalsuite_package_AUTHORS_deps])
AC_OUTPUT
AC_MSG_NOTICE([
==================================================
LALSuite has now been successfully configured:

* LALFrame library support is $LALFRAME_ENABLE_VAL
* LALMetaIO library support is $LALMETAIO_ENABLE_VAL
* LALSimulation library support is $LALSIMULATION_ENABLE_VAL
* LALBurst library support is $LALBURST_ENABLE_VAL
* LALInspiral library support is $LALINSPIRAL_ENABLE_VAL
* LALInference library support is $LALINFERENCE_ENABLE_VAL
* LALPulsar library support is $LALPULSAR_ENABLE_VAL
* LALApps library support is $LALAPPS_ENABLE_VAL
* Doxygen documentation is $DOXYGEN_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALSuite,
and run 'make install' to install LALSuite.
==================================================
])
