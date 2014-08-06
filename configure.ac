AC_PREREQ([2.63])
AC_INIT([LALSuite],[6.23.1],[lal-discuss@ligo.org])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_PREFIX_DEFAULT(["`pwd`/_inst"])
AC_CONFIG_FILES([Makefile])
AM_INIT_AUTOMAKE([1.11 foreign subdir-objects color-tests parallel-tests])

LALSUITE_DISTCHECK_CONFIGURE_FLAGS

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

# directories of top-level build and source
lalsuite_abs_top_builddir="`pwd`"
#lalsuite_top_builddir='$(top_builddir)/..'
lalsuite_top_srcdir='$(top_srcdir)/..'

# configure a lalsuite package
AC_DEFUN([lalsuite_config_subdir],[
  m4_pushdef([lowercase],m4_translit([[$1]], [A-Z], [a-z]))
  m4_pushdef([uppercase],m4_translit([[$1]], [a-z], [A-Z]))

  # if $1 is enabled
  AS_IF([test "x${lowercase}" = xtrue],[

    # export precious environment variables
    uppercase[]_LIBS="${lalsuite_abs_top_builddir}/lowercase/src/lib[]lowercase.la"
    uppercase[]_CFLAGS="-I${lalsuite_abs_top_builddir}/lowercase/include"
    uppercase[]_DATA_PATH="${lalsuite_top_srcdir}/lowercase/test"
    uppercase[]_OCTAVE_PATH="${lalsuite_abs_top_builddir}/lowercase/octave"
    uppercase[]_PYTHON_PATH="${lalsuite_abs_top_builddir}/lowercase/python"
    export uppercase[]_LIBS uppercase[]_CFLAGS
    export uppercase[]_DATA_PATH uppercase[]_OCTAVE_PATH uppercase[]_PYTHON_PATH

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

# doxygen, lal, and lalsupport are always configured
AC_CONFIG_SUBDIRS(doxygen)
AC_CONFIG_SUBDIRS(lal)
LAL_LIBS="${lalsuite_abs_top_builddir}/lal/lib/lal/liblal.la"
LAL_CFLAGS="-I${lalsuite_abs_top_builddir}/lal/include"
LAL_DATA_PATH="${lalsuite_top_srcdir}/lal/test"
LAL_OCTAVE_PATH="${lalsuite_abs_top_builddir}/lal/octave"
LAL_PYTHON_PATH="${lalsuite_abs_top_builddir}/lal/python"
LALSUPPORT_LIBS="${lalsuite_abs_top_builddir}/lal/lib/lalsupport/src/liblalsupport.la"
LALSUPPORT_CFLAGS="-I${lalsuite_abs_top_builddir}/lal/include"
LALSUPPORT_DATA_PATH="${LAL_DATA_PATH}"
LALSUPPORT_OCTAVE_PATH="${LAL_OCTAVE_PATH}"
LALSUPPORT_PYTHON_PATH="${LAL_PYTHON_PATH}"
LALSUITE_BUILD="true"
export LALSUITE_BUILD
export LAL_LIBS LAL_CFLAGS LAL_DATA_PATH LAL_OCTAVE_PATH LAL_PYTHON_PATH
export LALSUPPORT_LIBS LALSUPPORT_CFLAGS LALSUPPORT_DATA_PATH LALSUPPORT_OCTAVE_PATH LALSUPPORT_PYTHON_PATH

# configure optional packages
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

# optionally configure lalapps
AS_IF([test "x$lalapps" = xtrue],[
  AC_CONFIG_SUBDIRS(lalapps)
  LALAPPS_ENABLE_VAL=ENABLED
],[
  LALAPPS_ENABLE_VAL=DISABLED
])

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
