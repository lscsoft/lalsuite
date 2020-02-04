#!/bin/bash

set -ex

LALINFERENCE_DIR=$(cd "$(dirname $(readlink -f "${BASH_SOURCE[0]}"))"/.. && pwd)
export USER=albert.einstein
export OMP_NUM_THREADS=1
sed \
  -e 's|/cvmfs/oasis.opensciencegrid.org/ligo/sw/conda/envs/igwn-py37/| |'\
  -e "/ligo-skymap-from-samples=/c\ligo-skymap-from-samples=/bin/true" \
  -e "/ligo-skymap-plot=/c\ligo-skymap-plot=/bin/true" \
  -e "/email=/c\email="test" " \
  -e 's|albert.einstein/public_html/|'$USER'/public_html/LVC/|' \
  -e 's|^#.*fake-cache=|fake-cache=|' \
  -e 's|^#.*dataseed|dataseed|' \
  -e 's|^#.*randomseed|randomseed|' \
  -e 's|^#.*0noise=|0noise=|' \
  -e 's|^#.*notification=Complete|notification=Complete|' \
  -e "/approx=SEOBNRv2_ROM_DoubleSpinpseudoFourPN,IMRPhenomPv2pseudoFourPN/c\approx=IMRPhenomPv2pseudoFourPN" \
  -e 's|roq=False|roq=True|' \
  -e 's|/home/cbc|'$PWD'|' \
  -e "/nparallel=/c\nparallel=1" \
  -e "/tolerance=/c\tolerance=1.0" \
  -e 's|nlive=512|nlive=256\
adapt-tau=3\
distance-min=1000\
distance-max=1500|' \
  -e "/mpirun=/c\mpirun=/bin/true " \
  -e 's|mpi_task_count=8|mpi_task_count=4|' \
  -e 's|machine-count=8|machine-count=4|' \
  -e 's|ntemps=8|ntemps=4|' \
  -e "/accounting_group=/c\accounting_group=ligo.dev.o3.cbc.pe.lalinference" \
  -e 's|sharedfs=False|sharedfs=True|' \
  -e 's|^resume=|#resume=|' \
  -e 's|^#.*mtotal-min=2|mtotal-min=55|' \
  -e 's|^#.*mtotal-max=35|mtotal-max=65|' \
  -e '/^#.*maxmcmc=/c\maxmcmc=100' \
  ${LALINFERENCE_DIR}/lib/lalinference_pipe_example.ini > example.ini
lalinference_pipe --run-path ./example -I lalinference/test/injection_standard.xml --daglog-path ./daglog ./example.ini
cd example/4s
time bash -ex ./lalinference_441417463-441417627.sh
