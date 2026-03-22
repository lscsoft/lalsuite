# syntax=docker/dockerfile:1

ARG UBU_VERSION

FROM ubuntu:${UBU_VERSION}

ARG UBU_VERSION

LABEL name="LALSuite CI Image - Ubuntu ${UBU_VERSION}"
LABEL maintainer="LALSuite Maintainers <lal-discuss@ligo.org>"
LABEL support="Best Effort"

SHELL ["/bin/bash", "-c"]

# run debconf noninteractively
ENV DEBIAN_FRONTEND=noninteractive

RUN <<EOF
set -ex

# update APT cache
apt-get -y -q update

# upgrade distribution
apt-get -y -q upgrade

# install required packages
apt-get -y -q install \
    apt-transport-https \
    apt-utils \
    build-essential \
    ccache \
    curl \
    devscripts \
    dpkg-dev \
    equivs \
    git \
    git-lfs \
    gpg \
    lintian \
    local-apt-repository \
    lsb-release \
    python3-pip \
    xz-utils \
    ;
python3 -m pip install --upgrade --force-reinstall --ignore-installed --break-system-packages pip
python3 -m pip install --break-system-packages 'lintian-codeclimate>=0.1.1'

# install Git LFS
git lfs install

# cleanup
apt-get clean

EOF

# add lscsoft repository
COPY <<EOF /etc/apt/sources.list.d/lscsoft.list
deb [signed-by=/etc/apt/trusted.gpg.d/hypatia-key.asc] https://hypatia.aei.mpg.de/ubuntu/ ${UBU_VERSION} lvk
EOF

RUN <<EOF
set -ex

# set up lscsoft repository
cat /etc/apt/sources.list.d/lscsoft.list
curl -O https://hypatia.aei.mpg.de/ubuntu/keys/hypatia-keyring.deb
dpkg -i hypatia-keyring.deb
rm -f hypatia-keyring.deb

# update APT cache
apt-get -y -q update

# upgrade distribution
apt-get -y -q upgrade

# print info
dpkg-query --list
python3 -m pip list

# cleanup
apt-get clean
python3 -m pip cache purge

EOF
