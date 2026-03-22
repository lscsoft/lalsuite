# ----------------------------------------------------------------------
# LALSuite GitLab-CI: initialise build for Conda packages
# ----------------------------------------------------------------------

# initialise Conda
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate base

# recreate environment for CI jobs
# - if not running under the Conda container
# - assume base environment is correctly set up
if ! conda list --name lalsuite-ci; then

    # configure Conda options
    conda config --set always_yes yes
    conda config --add channels conda-forge
    conda config --set channel_priority strict

    # create environment for CI jobs
    conda create --quiet --name lalsuite-ci

fi

# activate environment
conda activate lalsuite-ci

# print info
conda info --all
conda config --show-sources
conda config --show
conda list --name base
conda list --name lalsuite-ci
