HERE=$( cd "$(dirname $(readlink -f "${BASH_SOURCE[0]}" ) )"/.. && pwd)
export USER=albert.einstein
export OMP_NUM_THREADS=1
sed \
  -e 's|/home/albert.einstein/opt/lalsuite/| |'\
  -e "/ligo-skymap-from-samples=/c\ligo-skymap-from-samples=/bin/true" \
  -e "/ligo-skymap-plot=/c\ligo-skymap-plot=/bin/true" \
  -e "/email=/c\email="test" " \
  -e 's|albert.einstein/public_html/|'$USER'/public_html/LVC/|' \
  -e 's|#fake-cache|fake-cache|' \
  -e 's|# dataseed|dataseed|' \
  -e 's|#0noise=|0noise=|' \
  -e "/approx=SEOBNRv2_ROM_DoubleSpinpseudoFourPN,IMRPhenomPv2pseudoFourPN/c\approx=IMRPhenomPv2pseudoFourPN" \
  -e 's|roq=False|roq=True|' \
  -e 's|/home/cbc|'$PWD'|' \
  -e "/nparallel=/c\nparallel=1" \
  -e 's|tolerance=0.1|tolerance=100|' \
  -e 's|nlive=512|nlive=256\
adapt-tau=3|' \
  -e 's|# randomseed=4321|randomseed=4321|' \
  -e "/mpirun=/c\mpirun=/bin/true " \
  -e 's|mpi_task_count=8|mpi_task_count=4|' \
  -e 's|machine-count=8|machine-count=4|' \
  -e 's|ntemps=8|ntemps=4|' \
  -e 's|#notification=Complete|notification=Complete|' \
  -e "/accounting_group=/c\accounting_group=ligo.dev.o3.cbc.pe.lalinference" \
  ${HERE}/src/lalinference_pipe_example.ini > example.ini
lalinference_pipe --run-path ./example -I lalinference/test/injection_standard.xml --daglog-path ./daglog ./example.ini
cd example/4s
time ./lalinference_441417463-441417627.sh
