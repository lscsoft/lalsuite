#!/bin/sh -x

## Perform the pre-install for LAL.  The following commands will download,
## build, and install the software required to build LAL.  The software
## includes the following:
##verse
##	pkgconfig-0.15.0
##	fftw-3.0.1 (both single- and double-precision libraries required)
##	gsl-1.5
##	libframe-6.14 (optional but recommended)
##	libmetaio-5.4 (optional but recommended)
##/verse
## If this software is already on your system, you can use the existing
## software.  If some of the software is missing, you can use the appropriate
## part of these instructions to install that component.  Pre-compiled versions
## of the software are also available for installation.  See the RPMs
## section for instructions on obtaining these.
##
## The software is installed in the directory "LSCSOFT_LOCATION".
## If this variable is not set, it will be installed in "$HOME/opt/lscsoft"
## by default. To install in some other location, set "LSCSOFT_LOCATION"
## to that location.
## 
## The commands listed below are appropriate for a Bourne-shell (e.g., bash);
## they will need to be modified appropriately for C-shells (e.g., tcsh).

#verbatim
LSCSOFT_PREFIX=${LSCSOFT_LOCATION:-"$HOME/opt/lscsoft"}
LSCSOFT_INCDIR=$LSCSOFT_PREFIX/include
LSCSOFT_LIBDIR=$LSCSOFT_PREFIX/lib
LSCSOFT_SRCDIR=$LSCSOFT_PREFIX/src
LSCSOFT_ETCDIR=$LSCSOFT_PREFIX/etc
LSCSOFT_TMPDIR=$LSCSOFT_PREFIX/tmp
#/verbatim

## This is where to get sources:
#verbatim
LALSRCURL=http://www.lsc-group.phys.uwm.edu/lal/sources
#/verbatim

### uncomment to use lynx instead of curl
###curl() {
###lynx -dump $1
###}

### uncomment to use wget instead of curl
###curl() {
###wget -O- $1
###}

###
### the rest of this script should not need to be edited
###

### simple failure
#ignore
fail() {
  echo "!!! Failure" 1>&2
  exit 1
}
#/ignore

# setup directories
#verbatim
mkdir -p $LSCSOFT_PREFIX || fail
mkdir -p $LSCSOFT_INCDIR || fail
mkdir -p $LSCSOFT_LIBDIR || fail
mkdir -p $LSCSOFT_SRCDIR || fail
mkdir -p $LSCSOFT_ETCDIR || fail
mkdir -p $LSCSOFT_TMPDIR || fail
#/verbatim

# get required autoconf, automake, fftw3, frame, gsl, and metaio
# you can use "lynx -dump" or "wget -O-" instead of "curl"
#verbatim
curl $LALSRCURL/pkgconfig-0.15.0.tar.gz > $LSCSOFT_TMPDIR/pkgconfig-0.15.0.tar.gz || fail
curl $LALSRCURL/fftw-3.0.1.tar.gz > $LSCSOFT_TMPDIR/fftw-3.0.1.tar.gz || fail
curl $LALSRCURL/gsl-1.5.tar.gz > $LSCSOFT_TMPDIR/gsl-1.5.tar.gz || fail
curl $LALSRCURL/libframe-6.14.tar.gz > $LSCSOFT_TMPDIR/libframe-6.14.tar.gz || fail
curl $LALSRCURL/libmetaio-5.4.tar.gz > $LSCSOFT_TMPDIR/libmetaio-5.4.tar.gz || fail
#/verbatim

# unpack these archives in "LSCSOFT_SRCDIR"
#verbatim
cd $LSCSOFT_SRCDIR || fail
tar -zxvf $LSCSOFT_TMPDIR/pkgconfig-0.15.0.tar.gz || fail
tar -zxvf $LSCSOFT_TMPDIR/fftw-3.0.1.tar.gz || fail
tar -zxvf $LSCSOFT_TMPDIR/gsl-1.5.tar.gz || fail
tar -zxvf $LSCSOFT_TMPDIR/libframe-6.14.tar.gz || fail
tar -zxvf $LSCSOFT_TMPDIR/libmetaio-5.4.tar.gz || fail
#/verbatim

# build and install pkg-config
#verbatim
cd $LSCSOFT_SRCDIR/pkgconfig-0.15.0 || fail
./configure --prefix=$LSCSOFT_PREFIX || fail
make || fail
make install || fail
#/verbatim

# build and install fftw3
#verbatim
cd $LSCSOFT_SRCDIR/fftw-3.0.1 || fail
./configure --prefix=$LSCSOFT_PREFIX --enable-shared --enable-float || fail
make  # note: ignore fail... the build fails on MacOSX, but not seriously
make install # note: ignore fail
make distclean || fail
./configure --prefix=$LSCSOFT_PREFIX --enable-shared || fail
make # note: ignore fail
make install # note: ignore fail
#/verbatim

# build and install gsl
#verbatim
cd $LSCSOFT_SRCDIR/gsl-1.5 || fail
./configure --prefix=$LSCSOFT_PREFIX || fail
make || fail
make install || fail
#/verbatim

# build and install libframe
#verbatim
cd $LSCSOFT_SRCDIR/libframe-6.14 || fail
./configure --prefix=$LSCSOFT_PREFIX --disable-octave || fail
make || fail
make install || fail
#/verbatim

# build and install libmetaio
#verbatim
cd $LSCSOFT_SRCDIR/libmetaio-5.4 || fail
./configure --prefix=$LSCSOFT_PREFIX || fail
make || fail
make install || fail
#/verbatim

### write environment configuration file
#ignore
rm -f $LSCSOFT_ETCDIR/lscsoft-user-env.sh || fail
cat > $LSCSOFT_ETCDIR/lscsoft-user-env.sh <<\EOF
# Source this file to set up your environment to use lscsoft software.
# This requires that LSCSOFT_LOCATION be set.
# LSCSOFT_PREFIX will be set by this script to save the current location
# so that the old LSCSOFT_PREFIX information can be removed from your
# environment if LSCSOFT_LOCATION is changed and this file is resourced.
# If LSCSOFT_LOCATION is set but empty then the previous location is
# removed from the environment.

if [ "${LSCSOFT_LOCATION-X}" = "X" ]; then
  echo "ERROR: environment variable LSCSOFT_LOCATION not defined" 1>&2
  return 1
fi

if [ -n "${LSCSOFT_PREFIX}" ]; then
  PATH=`echo "${PATH}" | sed -e "s%:${LSCSOFT_PREFIX}[^:]*%%g" -e "s%^${LSCSOFT_PREFIX}[^:]*:\{0,1\}%%"`
  LD_LIBRARY_PATH=`echo "${LD_LIBRARY_PATH}" | sed -e "s%:${LSCSOFT_PREFIX}[^:]*%%g" -e "s%^${LSCSOFT_PREFIX}[^:]*:\{0,1\}%%"`
  MANPATH=`echo "${MANPATH}" | sed -e "s%:${LSCSOFT_PREFIX}[^:]*%%g" -e "s%^${LSCSOFT_PREFIX}[^:]*:\{0,1\}%%"`
  export PATH MANPATH LD_LIBRARY_PATH
fi

LSCSOFT_PREFIX=${LSCSOFT_LOCATION}
export LSCSOFT_PREFIX

if [ -n "${LSCSOFT_PREFIX}" ]; then
  PATH=`echo "${PATH}" | sed -e "s%:${LSCSOFT_PREFIX}[^:]*%%g" -e "s%^${LSCSOFT_PREFIX}[^:]*:\{0,1\}%%"`
  LD_LIBRARY_PATH=`echo "${LD_LIBRARY_PATH}" | sed -e "s%:${LSCSOFT_PREFIX}[^:]*%%g" -e "s%^${LSCSOFT_PREFIX}[^:]*:\{0,1\}%%"`
  MANPATH=`echo "${MANPATH}" | sed -e "s%:${LSCSOFT_PREFIX}[^:]*%%g" -e "s%^${LSCSOFT_PREFIX}[^:]*:\{0,1\}%%"`
  PATH=${LSCSOFT_LOCATION}/bin:${PATH}
  LD_LIBRARY_PATH=${LSCSOFT_LOCATION}/lib:${LD_LIBRARY_PATH}
  MANPATH=${LSCSOFT_LOCATION}/man:${MANPATH}
  export PATH MANPATH LD_LIBRARY_PATH
fi
EOF
#/ignore
#ignore
rm -f $LSCSOFT_ETCDIR/lscsoft-user-env.csh || fail
cat > $LSCSOFT_ETCDIR/lscsoft-user-env.csh <<\EOF
# Source this file to set up your environment to use lscsoft software.
# This requires that LSCSOFT_LOCATION be set.
# LSCSOFT_PREFIX will be set by this script to save the current location
# so that the old LSCSOFT_PREFIX information can be removed from your
# environment if LSCSOFT_LOCATION is changed and this file is resourced.
# If LSCSOFT_LOCATION is set but empty then the previous location is
# removed from the environment.

if ( ! ${?LSCSOFT_LOCATION} ) then
  echo "ERROR: environment variable LSCSOFT_LOCATION not defined"
  exit 1
endif

if ( ! ${?LD_LIBRARY_PATH} ) then
  setenv LD_LIBRARY_PATH ''
endif

if ( ! ${?MANPATH} ) then
  setenv MANPATH ''
endif


if ( ${?LSCSOFT_PREFIX} ) then
  if (  "${LSCSOFT_PREFIX}" != "" ) then
    setenv PATH `echo "${PATH}" | sed -e "s%:${LSCSOFT_PREFIX}[^:]*%%g" -e "s%^${LSCSOFT_PREFIX}[^:]*:\{0,1\}%%"`
    setenv LD_LIBRARY_PATH `echo "${LD_LIBRARY_PATH}" | sed -e "s%:${LSCSOFT_PREFIX}[^:]*%%g" -e "s%^${LSCSOFT_PREFIX}[^:]*:\{0,1\}%%"`
    setenv MANPATH `echo "${MANPATH}" | sed -e "s%:${LSCSOFT_PREFIX}[^:]*%%g" -e "s%^${LSCSOFT_PREFIX}[^:]*:\{0,1\}%%"`
  endif
endif

setenv LSCSOFT_PREFIX ${LSCSOFT_LOCATION}

if ( "${LSCSOFT_PREFIX}" != "" ) then
  setenv PATH `echo "${PATH}" | sed -e "s%:${LSCSOFT_PREFIX}[^:]*%%g" -e "s%^${LSCSOFT_PREFIX}[^:]*:\{0,1\}%%"`
  setenv LD_LIBRARY_PATH `echo "${LD_LIBRARY_PATH}" | sed -e "s%:${LSCSOFT_PREFIX}[^:]*%%g" -e "s%^${LSCSOFT_PREFIX}[^:]*:\{0,1\}%%"`
  setenv MANPATH `echo "${MANPATH}" | sed -e "s%:${LSCSOFT_PREFIX}[^:]*%%g" -e "s%^${LSCSOFT_PREFIX}[^:]*:\{0,1\}%%"`
  setenv PATH ${LSCSOFT_LOCATION}/bin:${PATH}
  setenv LD_LIBRARY_PATH ${LSCSOFT_LOCATION}/lib:${LD_LIBRARY_PATH}
  setenv MANPATH ${LSCSOFT_LOCATION}/man:${MANPATH}
endif
EOF
#/ignore

## To setup your environment to use the software that has been installed
## please add the following to your .profile if you use a bourne shell
## (e.g. bash):
### (don't actually do it... for illustration purposes only)
#ignore
if 0 ; then
#/ignore
#verbatim
LSCSOFT_LOCATION=${HOME}/opt/lscsoft # <- change this as appropriate
export LSCSOFT_LOCATION
if [ -f ${LSCSOFT_LOCATION}/etc/lscsoft-user-env.sh ] ; then
  . ${LSCSOFT_LOCATION}/etc/lscsoft-user-env.sh
fi
#ignore
fi
if 0 ; then
#/ignore
#/verbatim
## If you are using a C shell (e.g., tcsh), instead add these lines to
## your .login:
#verbatim
setenv LSCSOFT_LOCATION ${HOME}/opt/lscsoft # <- change this as appropriate
if ( -r ${LSCSOFT_LOCATION}/etc/lscsoft-user-env.csh ) then
  source ${LSCSOFT_LOCATION}/etc/lscsoft-user-env.csh
endif
#/verbatim
#ignore
fi
fi
#/ignore

### Print a message alerting user to set "LSCSOFT_LOCATION"
#ignore
cat <<EOF
=======================================================================

To setup your environment to use the software that has been installed
please add the following to your .profile:

  LSCSOFT_LOCATION=$LSCSOFT_LOCATION
  export LSCSOFT_LOCATION
  if [ -f \${LSCSOFT_LOCATION}/etc/lscsoft-user-env.sh ] ; then
    . \${LSCSOFT_LOCATION}/etc/lscsoft-user-env.sh
  fi

If you are using a C shell (e.g., tcsh), instead add these lines to
your .login:

  setenv LSCSOFT_LOCATION $LSCSOFT_LOCATION
  if ( -r \${LSCSOFT_LOCATION}/etc/lscsoft-user-env.csh ) then
    source \${LSCSOFT_LOCATION}/etc/lscsoft-user-env.csh
  endif

=======================================================================
EOF
#/ignore

### all done
#ignore
exit 0
#/ignore
