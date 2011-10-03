# Source this file to access LALApps
setenv LALAPPS_PREFIX /opt/lscsoft/lalapps
#setenv LALAPPS_PREFIX /usr
setenv PATH "${LALAPPS_PREFIX}/bin:${PATH}"
setenv PYSITE_PATH python`python -V |& cut -d' ' -f2 | cut -d. -f-2`
if ( $?MANPATH ) then
  setenv MANPATH "${LALAPPS_PREFIX}/share/man:${MANPATH}"
else
  setenv MANPATH "${LALAPPS_PREFIX}/share/man"
endif
if ( $?PYTHONPATH ) then
  setenv PYTHONPATH "${LALAPPS_PREFIX}/lib/${PYSITE_PATH}/site-packages:${PYTHONPATH}"
else
  setenv PYTHONPATH "${LALAPPS_PREFIX}/lib/${PYSITE_PATH}/site-packages"
endif
