# ----------------------------------------------------------------------
# LALSuite GitLab-CI: install Conda dependencies for top-level build
# ----------------------------------------------------------------------

# install LALSuite build dependencies
conda env update --quiet --name lalsuite-ci --file "${CI_PROJECT_DIR}/conda-dev-env.yml"

# activate environment
conda activate lalsuite-ci

# print info
conda list --name lalsuite-ci
