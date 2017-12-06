AC_PREREQ([2.63])
AC_INIT([LALSuite],[6.38.1],[lal-discuss@ligo.org])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_PREFIX_DEFAULT(["`pwd`/_inst"])
AC_CONFIG_FILES([Makefile])
AM_INIT_AUTOMAKE([ \
  1.11 \
  foreign \
  subdir-objects \
  color-tests \
  parallel-tests \
  tar-ustar \
  dist-xz \
  no-dist-gzip \
])

LALSUITE_DISTCHECK_CONFIGURE_FLAGS

# use silent build rules
AM_SILENT_RULES([yes])

# provide LAL library enable/disable options
LALSUITE_ENABLE_ALL_LAL
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

# configure Doxygen top-level page
LALSUITE_ENABLE_DOXYGEN
LALSUITE_ENABLE_MODULE([DOXYGEN])
AS_IF([test "x${doxygen}" = xtrue],[
  AC_PROG_SED
  AC_CONFIG_FILES([index.html:gnuscripts/lalsuite_index.html.in])
  AC_SUBST([DOXYGEN_PACKAGE_LINKS],[])
  AM_SUBST_NOTMAKE([DOXYGEN_PACKAGE_LINKS])
  AC_SUBST([DOXYGEN_LAST_PACKAGE],[])
  AM_SUBST_NOTMAKE([DOXYGEN_LAST_PACKAGE])
])
AC_DEFUN([lalsuite_config_doxygen],[
  m4_pushdef([lowercase],m4_translit([[$1]], [A-Z], [a-z]))
  AS_IF([test "x${doxygen}" = xtrue],[
    DOXYGEN_PACKAGE_LINKS="${DOXYGEN_PACKAGE_LINKS}"' <a target="docframe" href="lowercase/doxygen/out/index.html">$1</a>'
    DOXYGEN_LAST_PACKAGE="lowercase/doxygen/out/index.html"
  ])
  m4_popdef([lowercase])
])

# directories of top-level build and source
lalsuite_abs_top_builddir="`pwd`"
lalsuite_abs_top_srcdir='$(abs_top_srcdir)/..'

# configure a lalsuite package
AC_DEFUN([lalsuite_config_subdir],[
  m4_pushdef([lowercase],m4_translit([[$1]], [A-Z], [a-z]))
  m4_pushdef([uppercase],m4_translit([[$1]], [a-z], [A-Z]))

  # if $1 is enabled
  AS_IF([test "x${lowercase}" = xtrue],[

    # export precious environment variables
    uppercase[]_LIBS="${lalsuite_abs_top_builddir}/lowercase/src/lib[]lowercase.la"
    uppercase[]_CFLAGS="-I${lalsuite_abs_top_builddir}/lowercase/include"
    uppercase[]_DATA_PATH="${lalsuite_abs_top_srcdir}/lowercase/test"
    uppercase[]_OCTAVE_PATH="${lalsuite_abs_top_builddir}/lowercase/octave"
    uppercase[]_PYTHON_PATH="${lalsuite_abs_top_builddir}/lowercase/python"
    uppercase[]_HTMLDIR="${htmldir}/../lowercase"
    export uppercase[]_LIBS uppercase[]_CFLAGS uppercase[]_HTMLDIR
    export uppercase[]_DATA_PATH uppercase[]_OCTAVE_PATH uppercase[]_PYTHON_PATH

    # configure Doxygen top-level page
    lalsuite_config_doxygen([$1])

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
lalsuite_config_doxygen([LAL])
AC_CONFIG_SUBDIRS(lal)
LAL_LIBS="${lalsuite_abs_top_builddir}/lal/src/liblal.la"
LAL_CFLAGS="-I${lalsuite_abs_top_builddir}/lal/include"
LAL_DATA_PATH="${lalsuite_abs_top_srcdir}/lal/test"
LAL_OCTAVE_PATH="${lalsuite_abs_top_builddir}/lal/octave"
LAL_PYTHON_PATH="${lalsuite_abs_top_builddir}/lal/python"
LAL_HTMLDIR="${htmldir}/../lal"
LALSUPPORT_LIBS="${lalsuite_abs_top_builddir}/lal/src/support/liblalsupport.la"
LALSUPPORT_CFLAGS="-I${lalsuite_abs_top_builddir}/lal/include"
LALSUPPORT_DATA_PATH="${LAL_DATA_PATH}"
LALSUPPORT_OCTAVE_PATH="${LAL_OCTAVE_PATH}"
LALSUPPORT_PYTHON_PATH="${LAL_PYTHON_PATH}"
LALSUPPORT_HTMLDIR="${htmldir}/../lal"
LALSUITE_BUILD="true"
export LALSUITE_BUILD
export LAL_LIBS LAL_CFLAGS LAL_DATA_PATH LAL_OCTAVE_PATH LAL_PYTHON_PATH LAL_HTMLDIR
export LALSUPPORT_LIBS LALSUPPORT_CFLAGS LALSUPPORT_DATA_PATH LALSUPPORT_OCTAVE_PATH LALSUPPORT_PYTHON_PATH LALSUPPORT_HTMLDIR

# configure optional packages
lalsuite_config_subdir([LALFrame])
lalsuite_config_subdir([LALMetaIO])
lalsuite_config_subdir([LALXML])
lalsuite_config_subdir([LALSimulation])
lalsuite_config_subdir([LALBurst])
lalsuite_config_subdir([LALDetChar])
lalsuite_config_subdir([LALInspiral])
lalsuite_config_subdir([LALStochastic])
lalsuite_config_subdir([LALPulsar])
lalsuite_config_subdir([LALInference])

# optionally configure lalapps
AS_IF([test "x$lalapps" = xtrue],[
  lalsuite_config_doxygen([LALApps])
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
* Doxygen documentation is $DOXYGEN_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALSuite,
and run 'make install' to install LALSuite.
==================================================
])
