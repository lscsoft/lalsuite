# Source this file to set up your environment to use lscsoft software.
# This requires that LAL_LOCATION be set.
# LAL_PREFIX will be set by this script to save the current location
# so that the old LAL_PREFIX information can be removed from your
# environment if LAL_LOCATION is changed and this file is resourced.
# If LAL_LOCATION is set but empty then the previous location is
# removed from the environment.

if [ "${LAL_LOCATION-X}" = "X" ]; then
  echo "ERROR: environment variable LAL_LOCATION not defined" 1>&2
  return 1
fi

if [ -n "${LAL_PREFIX}" ]; then
  PATH=`echo "${PATH}" | sed -e "s%:${LAL_PREFIX}[^:]*%%g" -e "s%^${LAL_PREFIX}[^:]*:\{0,1\}%%"`
  LD_LIBRARY_PATH=`echo "${LD_LIBRARY_PATH}" | sed -e "s%:${LAL_PREFIX}[^:]*%%g" -e "s%^${LAL_PREFIX}[^:]*:\{0,1\}%%"`
  MANPATH=`echo "${MANPATH}" | sed -e "s%:${LAL_PREFIX}[^:]*%%g" -e "s%^${LAL_PREFIX}[^:]*:\{0,1\}%%"`
  PKG_CONFIG_PATH=`echo "${PKG_CONFIG_PATH}" | sed -e "s%:${LAL_PREFIX}[^:]*%%g" -e "s%^${LAL_PREFIX}[^:]*:\{0,1\}%%"`
  export PATH MANPATH LD_LIBRARY_PATH PKG_CONFIG_PATH
fi

LAL_PREFIX=${LAL_LOCATION}
export LAL_PREFIX

if [ -n "${LAL_PREFIX}" ]; then
  PATH=`echo "${PATH}" | sed -e "s%:${LAL_PREFIX}[^:]*%%g" -e "s%^${LAL_PREFIX}[^:]*:\{0,1\}%%"`
  LD_LIBRARY_PATH=`echo "${LD_LIBRARY_PATH}" | sed -e "s%:${LAL_PREFIX}[^:]*%%g" -e "s%^${LAL_PREFIX}[^:]*:\{0,1\}%%"`
  MANPATH=`echo "${MANPATH}" | sed -e "s%:${LAL_PREFIX}[^:]*%%g" -e "s%^${LAL_PREFIX}[^:]*:\{0,1\}%%"`
  PKG_CONFIG_PATH=`echo "${PKG_CONFIG_PATH}" | sed -e "s%:${LAL_PREFIX}[^:]*%%g" -e "s%^${LAL_PREFIX}[^:]*:\{0,1\}%%"`
  PATH=${LAL_LOCATION}/bin:${PATH}
  LD_LIBRARY_PATH=${LAL_LOCATION}/lib:${LD_LIBRARY_PATH}
  MANPATH=${LAL_LOCATION}/man:${MANPATH}
  PKG_CONFIG_PATH=${LAL_LOCATION}/lib/pkgconfig:${PKG_CONFIG_PATH}
  export PATH MANPATH LD_LIBRARY_PATH PKG_CONFIG_PATH
fi
