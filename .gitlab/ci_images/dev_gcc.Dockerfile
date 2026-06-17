# syntax=docker/dockerfile:1

ARG GCC_VERSION
ARG DEB_VERSION

FROM gcc:${GCC_VERSION}-${DEB_VERSION}

ARG GCC_VERSION
ARG DEB_VERSION
ARG LCI_PKGLIST_X_LALAPPS
ARG TARBALL_NAME

LABEL name="LALSuite CI Image - GCC ${GCC_VERSION}"
LABEL maintainer="LALSuite Maintainers <lal-discuss@ligo.org>"
LABEL support="Best Effort"

SHELL ["/bin/bash", "-c"]

WORKDIR /tmp

# copy tarball for extraction of LALSuite build dependencies
COPY ./${TARBALL_NAME} .

# set compilers
ENV CC="/usr/local/bin/gcc"
ENV CXX="/usr/local/bin/g++"

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

RUN <<EOF
set -ex

# set up lscsoft repository
cat /etc/apt/sources.list.d/lscsoft.list
curl -O https://hypatia.aei.mpg.de/debian/keys/hypatia-keyring.deb
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

# print info
${CC} --version
${CXX} --version
dpkg-query --list
python3 -m pip list

# cleanup
apt-get clean
python3 -m pip cache purge
cd /tmp
rm -rf /tmp/*

EOF
