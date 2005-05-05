# Source this file to set up your environment to use lscsoft software.
# This requires that LAL_LOCATION be set.
# LAL_PREFIX will be set by this script to save the current location
# so that the old LAL_PREFIX information can be removed from your
# environment if LAL_LOCATION is changed and this file is resourced.
# If LAL_LOCATION is set but empty then the previous location is
# removed from the environment.

if ( ! ${?LAL_LOCATION} ) then
  echo "ERROR: environment variable LAL_LOCATION not defined"
  exit 1
endif

if ( ! ${?LD_LIBRARY_PATH} ) then
  setenv LD_LIBRARY_PATH ''
endif

if ( ! ${?MANPATH} ) then
  setenv MANPATH ''
endif

if ( ! ${?PKG_CONFIG_PATH} ) then
  setenv PKG_CONFIG_PATH ''
endif


if ( ${?LAL_PREFIX} ) then
  if (  "${LAL_PREFIX}" != "" ) then
    setenv PATH `echo "${PATH}" | sed -e "s%:${LAL_PREFIX}[^:]*%%g" -e "s%^${LAL_PREFIX}[^:]*:\{0,1\}%%"`
    setenv LD_LIBRARY_PATH `echo "${LD_LIBRARY_PATH}" | sed -e "s%:${LAL_PREFIX}[^:]*%%g" -e "s%^${LAL_PREFIX}[^:]*:\{0,1\}%%"`
    setenv MANPATH `echo "${MANPATH}" | sed -e "s%:${LAL_PREFIX}[^:]*%%g" -e "s%^${LAL_PREFIX}[^:]*:\{0,1\}%%"`
    setenv PKG_CONFIG_PATH `echo "${PKG_CONFIG_PATH}" | sed -e "s%:${LAL_PREFIX}[^:]*%%g" -e "s%^${LAL_PREFIX}[^:]*:\{0,1\}%%"`
  endif
endif

setenv LAL_PREFIX ${LAL_LOCATION}

if ( "${LAL_PREFIX}" != "" ) then
  setenv PATH `echo "${PATH}" | sed -e "s%:${LAL_PREFIX}[^:]*%%g" -e "s%^${LAL_PREFIX}[^:]*:\{0,1\}%%"`
  setenv LD_LIBRARY_PATH `echo "${LD_LIBRARY_PATH}" | sed -e "s%:${LAL_PREFIX}[^:]*%%g" -e "s%^${LAL_PREFIX}[^:]*:\{0,1\}%%"`
  setenv MANPATH `echo "${MANPATH}" | sed -e "s%:${LAL_PREFIX}[^:]*%%g" -e "s%^${LAL_PREFIX}[^:]*:\{0,1\}%%"`
  setenv PKG_CONFIG_PATH `echo "${PKG_CONFIG_PATH}" | sed -e "s%:${LAL_PREFIX}[^:]*%%g" -e "s%^${LAL_PREFIX}[^:]*:\{0,1\}%%"`
  setenv PATH ${LAL_LOCATION}/bin:${PATH}
  setenv LD_LIBRARY_PATH ${LAL_LOCATION}/lib:${LD_LIBRARY_PATH}
  setenv MANPATH ${LAL_LOCATION}/man:${MANPATH}
  setenv PKG_CONFIG_PATH ${LAL_LOCATION}/pkgconfig:${PKG_CONFIG_PATH}
endif
