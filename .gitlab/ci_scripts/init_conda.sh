# ----------------------------------------------------------------------
# LALSuite GitLab-CI: initialise build for Conda packages
# ----------------------------------------------------------------------

# initialise Conda
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate base

# configure Conda options
conda config --set always_yes yes
conda config --add channels conda-forge
conda config --set channel_priority strict
conda config --set remote_max_retries 5
conda config --set remote_backoff_factor 2
conda config --set remote_connect_timeout_secs 60.0
conda config --set remote_read_timeout_secs 180.0

if test -f /lalsuite-conda-container; then   # using a Conda container

    # activate environment
    conda activate lalsuite-ci

    # pin LALSuite build dependencies to prevent aggressive upgrades
    yq -r '.dependencies[] | select(type == "string")' "${CI_PROJECT_DIR}/conda-dev-env.yml" | while IFS= read -r pkg_dep; do
        pkg_name=$(echo "${pkg_dep}" | sed 's/[ <=>].*//')
        pkg_verspec=$(echo "${pkg_dep}" | sed "s/^${pkg_name} *//")
        pkg_version=$(conda list --name lalsuite-ci --full-name ${pkg_name} --json | yq -r '.[0].version')

        # if package has not maximum version restriction, use the currently installed version as a maximum
        case "${pkg_verspec}" in
            '='*)          # matches =, ==
                pkg_pin="${pkg_dep}"
                ;;
            *'<'*)         # matches <, <=
                pkg_pin="${pkg_dep}"
                ;;
            *'>'*)         # matches >, >=
                pkg_pin="${pkg_dep},<=${pkg_version}"
                ;;
            '')            # no version specifier
                pkg_pin="${pkg_name} <=${pkg_version}"
                ;;
            *)
                echo "ERROR: could not parse Conda version specifier '${pkg_verspec}'"
                exit 1
                ;;
        esac
        conda config --add pinned_packages "${pkg_pin}"
        echo "Added Conda package pin ${pkg_pin}"

    done

else   # using a restricted Conda environment, e.g. for macOS jobs

    # recreate environment for CI jobs
    # - assume base environment is correctly set up
    conda create --quiet --name lalsuite-ci
    conda activate lalsuite-ci

fi

# print info
conda info --all
conda config --show-sources
conda config --show
conda list --name base
conda list --name lalsuite-ci
