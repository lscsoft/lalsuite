# syntax=docker/dockerfile:1

FROM igwn/base:conda

ARG PYTHON_VERSION

LABEL name="LALSuite CI Image - Conda - Python ${PYTHON_VERSION}"
LABEL maintainer="LALSuite Maintainers <lal-discuss@ligo.org>"
LABEL support="Best Effort"

SHELL ["/bin/bash", "-c"]

# copy Conda environments
COPY ./conda-dev-env.yml .

# Conda root location
ENV CONDA_ROOT=/opt/conda

# Python version
ENV LCI_PYTHON_VERSION=${PYTHON_VERSION}

RUN <<EOF
set -ex

# initialise Conda
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate base

# configure Conda options
conda config --set always_yes yes
conda config --add channels conda-forge
conda config --set channel_priority strict

# update base environment
conda install --quiet --name base \
    "conda-build!=3.18.10,!=24.11.1" \
    "conda-smithy>=3.7.5" \
    "conda>=25" \
    conda-forge-pinning \
    conda-libmamba-solver \
    xz \
    ;

# create environment for CI jobs
conda create --quiet --name lalsuite-ci python=="${PYTHON_VERSION}"
conda activate lalsuite-ci

# pin Python version to PYTHON_VERSION
conda config --add pinned_packages python=="${PYTHON_VERSION}"

# install required CI packages
conda install --quiet --name lalsuite-ci \
    "python-gitlab>=5.6.0" \
    abi-compliance-checker \
    abi-dumper \
    chardet \
    charset-normalizer \
    coloredlogs \
    coverage \
    curl \
    gcovr \
    natsort \
    perl \
    pip \
    pipx \
    pre-commit \
    pycobertura \
    python-debian \
    pyyaml \
    requests \
    texlive-core \
    urllib3 \
    ;
python3 -m pip install \
    "py-api-dumper>=4.0" \
    pypi-cleanup \
    python-rpm-spec \
    twine \
    ;

# install LALSuite build dependencies
conda env update --quiet --name lalsuite-ci --file ./conda-dev-env.yml
rm -f ./conda-dev-env.yml

# print info
conda info --all
conda config --show-sources
conda config --show
conda list --name base
conda list --name lalsuite-ci

# cleanup
conda clean --all
python3 -m pip cache purge

EOF
