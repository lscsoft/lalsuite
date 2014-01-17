AC_PREREQ([2.63])
AC_INIT([LALSuite],[6.21.1],[lal-discuss@ligo.org])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_PREFIX_DEFAULT(["`pwd`/_inst"])
AC_CONFIG_FILES([Makefile])
AM_INIT_AUTOMAKE([1.11 foreign subdir-objects color-tests parallel-tests])

LALSUITE_DISTCHECK_CONFIGURE_FLAGS

# always prefer uninstalled LAL libraries
LALSUITE_PKG_SUFFIX="-uninstalled"
export LALSUITE_PKG_SUFFIX

# provide LAL library enable/disable options
LALSUITE_ENABLE_ALL_LAL
doxygen=true
lal=true
LALSUITE_ENABLE_LALFRAME
LALSUITE_ENABLE_LALMETAIO
LALSUITE_ENABLE_LALXML
LALSUITE_ENABLE_LALSIMULATION
LALSUITE_ENABLE_LALBURST
LALSUITE_ENABLE_LALDETCHAR
LALSUITE_ENABLE_LALINSPIRAL
LALSUITE_ENABLE_LALSTOCHASTIC
LALSUITE_ENABLE_LALPULSAR
LALSUITE_ENABLE_LALINFERENCE
LALSUITE_ENABLE_LALAPPS

# determine subdirectories to configure/build
AC_DEFUN([lalsuite_config_subdir],[
  m4_pushdef([lowercase],m4_translit([[$1]], [A-Z], [a-z]))
  m4_pushdef([uppercase],m4_translit([[$1]], [a-z], [A-Z]))

  # if $1 is enabled
  AS_IF([test "x${lowercase}" = xtrue],[

    # add $1 build directory to PKG_CONFIG_PATH
    LALSUITE_ADD_PATH(PKG_CONFIG_PATH,"`pwd`/lowercase")
    export PKG_CONFIG_PATH

    # configure $1
    AC_CONFIG_SUBDIRS(lowercase)

    # set enable/disable string
    uppercase[]_ENABLE_VAL=ENABLED
  ],[
    uppercase[]_ENABLE_VAL=DISABLED
  ])

  m4_popdef([lowercase])
  m4_popdef([uppercase])
])
lalsuite_config_subdir([doxygen])
lalsuite_config_subdir([lal])
lalsuite_config_subdir([lalframe])
lalsuite_config_subdir([lalmetaio])
lalsuite_config_subdir([lalxml])
lalsuite_config_subdir([lalsimulation])
lalsuite_config_subdir([lalburst])
lalsuite_config_subdir([laldetchar])
lalsuite_config_subdir([lalinspiral])
lalsuite_config_subdir([lalstochastic])
lalsuite_config_subdir([lalpulsar])
lalsuite_config_subdir([lalinference])
lalsuite_config_subdir([lalapps])

AC_SUBST([ac_configure_args])
AC_OUTPUT
AC_MSG_NOTICE([
==================================================
LALSuite has now been successfully configured:

* LALFrame library support is $LALFRAME_ENABLE_VAL
* LALMetaIO library support is $LALMETAIO_ENABLE_VAL
* LALXML library support is $LALXML_ENABLE_VAL
* LALSimulation library support is $LALSIMULATION_ENABLE_VAL
* LALBurst library support is $LALBURST_ENABLE_VAL
* LALDetChar library support is $LALDETCHAR_ENABLE_VAL
* LALInspiral library support is $LALINSPIRAL_ENABLE_VAL
* LALStochastic library support is $LALSTOCHASTIC_ENABLE_VAL
* LALPulsar library support is $LALPULSAR_ENABLE_VAL
* LALInference library support is $LALINFERENCE_ENABLE_VAL
* LALApps library support is $LALAPPS_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALSuite,
and run 'make install' to install LALSuite.
==================================================
])
