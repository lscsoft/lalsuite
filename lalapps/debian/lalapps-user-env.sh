# Source this file to access LALApps
PATH="/opt/lscsoft/lalapps/bin:${PATH}"
MANPATH="/opt/lscsoft/lalapps/share/man:${MANPATH}"
PYSITE_PATH=python`python -V 2>&1 | cut -d' ' -f2 | cut -d. -f-2`
PYTHONPATH="/opt/lscsoft/lalapps/lib/${PYSITE_PATH}/site-packages:${PYTHONPATH}"
export PATH MANPATH PYTHONPATH
