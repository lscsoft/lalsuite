# syntax=docker/dockerfile:1

ARG CUDA_VERSION
ARG UBU_VERSION
ARG UBU_VERSION_N

FROM nvidia/cuda:${CUDA_VERSION}-devel-ubuntu${UBU_VERSION_N}

ARG CUDA_VERSION
ARG UBU_VERSION
ARG UBU_VERSION_N
ARG LCI_PKGLIST_X_LALAPPS
ARG TARBALL_NAME

LABEL name="LALSuite CI Image - CUDA ${CUDA_VERSION} (Ubuntu ${UBU_VERSION})"
LABEL maintainer="LALSuite Maintainers <lal-discuss@ligo.org>"
LABEL support="Best Effort"

SHELL ["/bin/bash", "-c"]

WORKDIR /tmp

# copy tarball for extraction of LALSuite build dependencies
COPY ./${TARBALL_NAME} .

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

# update APT cache
apt-get -y -q update

# upgrade distribution
apt-get -y -q upgrade

# install LALSuite build dependencies
tar xf ./${TARBALL_NAME}
rm -f ./${TARBALL_NAME}
cd ./lalsuite-*
for subdir in ${LCI_PKGLIST_X_LALAPPS} lalapps; do
    pushd ${subdir}
    apt-get build-dep -y -q . || apt-get install -f || true
    popd
done

# install latest LALSuite release, if available
if ! apt-get -y -q install \
    $(printf "lib%s-dev " ${LCI_PKGLIST_X_LALAPPS}) \
    $(printf "python3-%s " ${LCI_PKGLIST_X_LALAPPS}) \
    ${LCI_PKGLIST_X_LALAPPS} lalapps
then
    echo 'WARNING: no LALSuite release available' > /no-lalsuite
    cat /no-lalsuite
fi

# print info
dpkg-query --list
python3 -m pip list

# cleanup
apt-get clean
python3 -m pip cache purge
cd /tmp
rm -rf /tmp/*

EOF
