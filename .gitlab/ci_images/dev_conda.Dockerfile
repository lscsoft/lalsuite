# syntax=docker/dockerfile:1

FROM igwn/base:conda

ARG PYTHON_VERSION
ARG LCI_PKGLIST_X_LALAPPS

LABEL name="LALSuite CI Image - Conda - Python ${PYTHON_VERSION}"
LABEL maintainer="LALSuite Maintainers <lal-discuss@ligo.org>"
LABEL support="Best Effort"

SHELL ["/bin/bash", "-c"]

WORKDIR /tmp

# copy Conda environments
COPY ./conda-dev-env.yml .

# Conda root location
ENV CONDA_ROOT=/opt/conda

# Python version
ENV LCI_PYTHON_VERSION=${PYTHON_VERSION}

# run debconf noninteractively
ENV DEBIAN_FRONTEND=noninteractive

# APT configuration
COPY <<EOF /etc/apt/apt.conf.d/99norecommends
APT::Install-Recommends "false";
APT::Install-Suggests "false";
EOF

RUN <<EOF
set -ex

# update APT cache
apt-get -y -q update

# upgrade distribution
apt-get -y -q upgrade

# install required packages
apt-get -y -q install \
    git \
    git-lfs \
    xz-utils \
    yq \
    ;

# install Git LFS
git lfs install

# cleanup
apt-get clean

EOF

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
    ;

# create environment for CI jobs with Python version pinned to ${PYTHON_VERSION}
conda create --quiet --name lalsuite-ci "python==${PYTHON_VERSION}"
conda config --add pinned_packages "python==${PYTHON_VERSION}"
conda activate lalsuite-ci

# install LALSuite build dependencies
conda env update --quiet --name lalsuite-ci --file ./conda-dev-env.yml

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

# create environment for testing package upgrade
conda create --quiet --name lalsuite-ci-upgrade

# install latest LALSuite release, if available
if ! conda install --quiet --name lalsuite-ci-upgrade \
    ${LCI_PKGLIST_X_LALAPPS} lalapps
then
    echo X > /no-lalsuite
    echo 'WARNING: no LALSuite release available, will skip package upgrade test'
fi

# create file to indicate a Conda container
# - as opposed to macOS jobs which use a restricted Conda environment
echo X > /lalsuite-conda-container

# print info
conda info --all
conda config --show-sources
conda config --show
conda list --name base
conda list --name lalsuite-ci
conda list --name lalsuite-ci-upgrade

# cleanup
conda clean --all
python3 -m pip cache purge
cd /tmp
rm -rf /tmp/*

EOF
