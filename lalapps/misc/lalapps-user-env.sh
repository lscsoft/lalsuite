# Source this file to set up your environment to use lscsoft software.
# This requires that LALAPPS_LOCATION be set.
# LALAPPS_PREFIX will be set by this script to save the current location
# so that the old LALAPPS_PREFIX information can be removed from your
# environment if LALAPPS_LOCATION is changed and this file is resourced.
# If LALAPPS_LOCATION is set but empty then the previous location is
# removed from the environment.

if [ "${LALAPPS_LOCATION-X}" = "X" ]; then
  echo "ERROR: environment variable LALAPPS_LOCATION not defined" 1>&2
  return 1
fi

if [ -n "${LALAPPS_PREFIX}" ]; then
  PATH=`echo "${PATH}" | sed -e "s%:${LALAPPS_PREFIX}[^:]*%%g" -e "s%^${LALAPPS_PREFIX}[^:]*:\{0,1\}%%"`
  MANPATH=`echo "${MANPATH}" | sed -e "s%:${LALAPPS_PREFIX}[^:]*%%g" -e "s%^${LALAPPS_PREFIX}[^:]*:\{0,1\}%%"`
  export PATH MANPATH 
fi

LALAPPS_PREFIX=${LALAPPS_LOCATION}
export LALAPPS_PREFIX

if [ -n "${LALAPPS_PREFIX}" ]; then
  PATH=`echo "${PATH}" | sed -e "s%:${LALAPPS_PREFIX}[^:]*%%g" -e "s%^${LALAPPS_PREFIX}[^:]*:\{0,1\}%%"`
  MANPATH=`echo "${MANPATH}" | sed -e "s%:${LALAPPS_PREFIX}[^:]*%%g" -e "s%^${LALAPPS_PREFIX}[^:]*:\{0,1\}%%"`
  PATH=${LALAPPS_LOCATION}/bin:${PATH}
  MANPATH=${LALAPPS_LOCATION}/man:${MANPATH}
  export PATH MANPATH 
fi
