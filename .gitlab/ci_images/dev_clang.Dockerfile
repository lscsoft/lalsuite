# syntax=docker/dockerfile:1

ARG DEB_VERSION

FROM debian:${DEB_VERSION}-slim

ARG CLANG_VERSION_NAME
ARG CLANG_VERSION_SUFFIX
ARG DEB_VERSION

LABEL name="LALSuite CI Image - Clang ${CLANG_VERSION_NAME}"
LABEL maintainer="LALSuite Maintainers <lal-discuss@ligo.org>"
LABEL support="Best Effort"

SHELL ["/bin/bash", "-c"]

# set compilers
ENV CC="clang${CLANG_VERSION_SUFFIX}"
ENV CXX="clang++${CLANG_VERSION_SUFFIX}"

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
    git \
    git-lfs \
    gpg \
    lsb-release \
    python3-pip \
    xz-utils \
    ;
python3 -m pip install --upgrade --force-reinstall --ignore-installed --break-system-packages pip

# install Git LFS
git lfs install

# cleanup
apt-get clean

EOF

# add lscsoft repository
COPY <<EOF /etc/apt/sources.list.d/lscsoft.list
deb [signed-by=/etc/apt/trusted.gpg.d/hypatia-key.asc] https://hypatia.aei.mpg.de/debian/ ${DEB_VERSION} lvk
EOF

# add llvm repository
COPY <<EOF /etc/apt/sources.list.d/llvm.list
deb [signed-by=/usr/share/keyrings/llvm-snapshot.gpg] http://apt.llvm.org/${DEB_VERSION}/ llvm-toolchain-${DEB_VERSION}${CLANG_VERSION_SUFFIX} main
EOF

RUN <<EOF
set -ex

# set up lscsoft repository
cat /etc/apt/sources.list.d/lscsoft.list
curl -O https://hypatia.aei.mpg.de/debian/keys/hypatia-keyring.deb
dpkg -i hypatia-keyring.deb
rm -f hypatia-keyring.deb

# set up llvm repository
curl -sSL https://apt.llvm.org/llvm-snapshot.gpg.key | gpg --dearmor --yes --output /usr/share/keyrings/llvm-snapshot.gpg

# update APT cache
apt-get -y -q update

# upgrade distribution
apt-get -y -q upgrade

# install Clang
apt-get -y -q install clang${CLANG_VERSION_SUFFIX}

# print info
dpkg-query --list
python3 -m pip list

# cleanup
apt-get clean
python3 -m pip cache purge

EOF
