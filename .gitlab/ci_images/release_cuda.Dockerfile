# syntax=docker/dockerfile:1

ARG CUDA_VERSION
ARG UBU_VERSION
ARG UBU_VERSION_N

FROM nvidia/cuda:${CUDA_VERSION}-runtime-ubuntu${UBU_VERSION_N}

ARG CUDA_VERSION
ARG UBU_VERSION
ARG UBU_VERSION_N

LABEL name="LALSuite Release Image - CUDA ${CUDA_VERSION} (Ubuntu ${UBU_VERSION})"
LABEL maintainer="LALSuite Maintainers <lal-discuss@ligo.org>"
LABEL support="Best Effort"

SHELL ["/bin/bash", "-c"]

# copy packages
COPY ./packages/ /packages/

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
    curl \
    gpg \
    local-apt-repository \
    ;

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

# install latest release
apt-get -y -q install lalsuite lalsuite-dev

# remove lalsuite meta-packages
dpkg -r lalsuite lalsuite-dev

# create local repo for upstream .debs
upstream_debs=$(find /packages -name '*.deb')
echo "===== upstream .debs"
printf "%s\n" ${upstream_debs}
echo "====="
mkdir -pv /srv/local-apt-repository
cp -v ${upstream_debs} /srv/local-apt-repository
/usr/lib/local-apt-repository/rebuild

# upgrade all packages
apt-get -y -q update
apt-get -y dist-upgrade

# print info
dpkg-query --list

# cleanup
apt-get clean
rm -rf /packages
rm -rf /srv/local-apt-repository

EOF
