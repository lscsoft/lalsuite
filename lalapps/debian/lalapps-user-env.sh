# Source this file to access LALApps
LALAPPS_PREFIX=/opt/lscsoft/lalapps
#LALAPPS_PREFIX=/usr
export LALAPPS_PREFIX
PATH="${LALAPPS_PREFIX}/bin:${PATH}"
MANPATH="${LALAPPS_PREFIX}/share/man:${MANPATH}"
PYSITE_PATH=python`python -V 2>&1 | cut -d' ' -f2 | cut -d. -f-2`
PYTHONPATH="${LALAPPS_PREFIX}/lib/${PYSITE_PATH}/site-packages:${PYTHONPATH}"
export PATH MANPATH PYTHONPATH
