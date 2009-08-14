# Source this file to access LALApps
setenv PATH "/opt/lscsoft/lalapps/bin:${PATH}"
setenv PYSITE_PATH python`python -V 2>&1 | cut -d' ' -f2 | cut -d. -f-2`
if ( $?MANPATH ) then
  setenv MANPATH "/opt/lscsoft/lalapps/share/man:${MANPATH}"
else
  setenv MANPATH "/opt/lscsoft/lalapps/share/man"
endif
if ( $?PYTHONPATH ) then
  setenv PYTHONPATH "/opt/lscsoft/lalapps/lib/${PYSITE_PATH}/site-packages:${PYTHONPATH}"
else
  setenv PYTHONPATH "/opt/lscsoft/lalapps/lib/${PYSITE_PATH}/site-packages"
endif
