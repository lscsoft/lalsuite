# ----------------------------------------------------------------------
# LALSuite GitLab-CI: make wheels
# ----------------------------------------------------------------------

# create virtual environments for building and testing
for venv in venv-build venv-test; do
    ${PYTHON} -m venv "${LCI_TMPDIR}/${venv}"
    source "${LCI_TMPDIR}/${venv}/bin/activate"
    ${LCI_SCRIPTS}/retry python -m pip install --upgrade pip
done

# activate building virtual environment
source "${LCI_TMPDIR}/venv-build/bin/activate"

# fall-back LAL data path
# - set ${base} to data path relative to binary LAL library install location
# - use the solar_system_ephemerides package for LALPulsar ephemeris files
sse="solar_system_ephemerides/ephemerides"
fallback_data_path='${base}/lalapps/data:${base}/${sse}/earth:${base}/${sse}/sun:${base}/${sse}/time'

os="${CI_JOB_GROUP_NAME##*:}"
echo "os=${os}"
case ${os} in

    # build recipe for standalone wheels on Linux
    linux-*)

        # build against numpy>=2.0.0rc1 for ABI compatibility with Numpy 1.x and 2.x.
        # - see https://numpy.org/devdocs/dev/depending_on_numpy.html#numpy-2-abi-handling.
        ${LCI_SCRIPTS}/retry python -m pip install \
            'numpy>=2.0.0' \
            setuptools \
            wheel
        python -m pip list installed

        # setup build environment
        source ${LCI_SCRIPTS}/build_env.sh

        # configure
        eval ./configure \
            PYTHON="${VIRTUAL_ENV}/bin/python" \
            --with-fallback-data-path="$(base=".."; eval echo ${fallback_data_path})" \
            --disable-doxygen \
            --enable-mpi \
            --without-ephem \
            "${LCI_CONFIGURE_FLAGS}"

        # build wheel
        make -j${CPU_COUNT} wheel

        # bundle and fix up dependent shared libraries
        auditwheel repair wheel/*.whl

        ;;

    # build recipe for standalone wheels on macOS
    macos-*)

        # build against numpy>=2.0.0rc1 for ABI compatibility with Numpy 1.x and 2.x.
        # - see https://numpy.org/devdocs/dev/depending_on_numpy.html#numpy-2-abi-handling.
        ${LCI_SCRIPTS}/retry python -m pip install \
            delocate \
            'numpy>=2.0.1' \
            setuptools \
            wheel
        python -m pip list installed

        # setup build environment
        source ${LCI_SCRIPTS}/build_env.sh

        # configure
        eval ./configure \
            PYTHON="${VIRTUAL_ENV}/bin/python" \
            LDFLAGS="-Wl,-headerpad_max_install_names" \
            --with-fallback-data-path="$(base="../.."; eval echo ${fallback_data_path})" \
            --disable-doxygen \
            --enable-mpi \
            --without-ephem \
            "${LCI_CONFIGURE_FLAGS}"

        # build wheel
        make -j${CPU_COUNT} wheel

        # bundle and fix up dependent shared libraries
        python -m delocate.cmd.delocate_wheel -v -w wheelhouse wheel/*.whl

        ;;

    *)
        echo "ERROR: unknown OS choice '${os}'"
        exit 1
        ;;

esac

# remove build tree first to make sure wheel does not
# access it, e.g. try to find LAL data files there
rm -rf wheel/build/

# activate testing virtual environment
source "${LCI_TMPDIR}/venv-test/bin/activate"

# install packages required for testing
${LCI_SCRIPTS}/retry python -m pip install \
    solar_system_ephemerides

# install wheel
python -m pip install wheelhouse/*

# check metadata
python -m pip show lalsuite

# verify wheels
# - skipped for now, see https://github.com/scipy/scipy/issues/21436
### python -m pip check lalsuite

# test loading Python modules and finding data files
env LAL_DEBUG_LEVEL=info python - <<EOF
import lal
import lalframe
import lalmetaio
import lalsimulation
import lalpulsar
series = lal.CreateREAL8FrequencySeries('psd', 0, 0, 1, None, 4096)
lalsimulation.SimNoisePSDaLIGOAPlusDesignSensitivityT1800042(series, 10)
lalpulsar.InitBarycenter("earth00-40-DE405.dat.gz", "sun00-40-DE405.dat.gz")
EOF

# test running Python scripts
lal_path2cache --help

# test running C executables and finding data files
lalapps_version
env LAL_DEBUG_LEVEL=info lalsim-detector-noise --official --aligo-nsnsopt --duration 1 >/dev/null
env LAL_DEBUG_LEVEL=info lalpulsar_PrintDetectorState --detector H1 --Alpha 4.65 --Delta -0.51 --timeGPS 1100000000 >/dev/null
