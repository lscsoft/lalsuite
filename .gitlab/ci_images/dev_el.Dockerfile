# syntax=docker/dockerfile:1

ARG EL_VERSION

FROM igwn/base:${EL_VERSION}

ARG EL_VERSION
ARG LCI_PKGLIST_X_LALAPPS
ARG TARBALL_NAME

LABEL name="LALSuite CI Image - Enterprise Linux ${EL_VERSION}"
LABEL maintainer="LALSuite Maintainers <lal-discuss@ligo.org>"
LABEL support="Reference Platform"

SHELL ["/bin/bash", "-c"]

WORKDIR /tmp

# copy tarball for extraction of LALSuite build dependencies
COPY ./${TARBALL_NAME} .

# write RPM macros
COPY <<EOF /root/.rpmmacros
%dist .${EL_VERSION}
EOF

RUN <<EOF
set -ex

# disable repos we don't use
DISABLE_REPOS="htcondor osg"
if [ "X${DISABLE_REPOS}" != X ]; then
    dnf -y -q install --disablerepo ${DISABLE_REPOS/ /,} "dnf-command(config-manager)"
    dnf config-manager --disable ${DISABLE_REPOS}
else
    dnf -y -q install "dnf-command(config-manager)"
fi

# upgrade distribution
dnf -y -q upgrade

# install development tools
dnf -y -q group install "Development Tools"

# install required packages
dnf -y -q install \
    createrepo \
    curl \
    findutils \
    igwn-packaging-tools \
    rpm-build \
    rpmlint \
    xz \
    ;
python3 -m pip install --upgrade pip
python3 -m pip install rpmlint-codeclimate

# install LALSuite build dependencies
tar xf ./${TARBALL_NAME}
rm -f ./${TARBALL_NAME}
cd ./lalsuite-*
specfiles=
for subdir in ${LCI_PKGLIST_X_LALAPPS} lalapps; do
    specfile=$(ls -1 ${subdir}/*.spec)
    specfiles="${specfiles} ${specfile}"
done
dnf builddep -y -q ${specfiles} || true

# install latest LALSuite release, if available
if ! dnf -y -q install \
    $(printf "lib%s-devel " ${LCI_PKGLIST_X_LALAPPS}) \
    $(printf "python3-%s " ${LCI_PKGLIST_X_LALAPPS}) \
    ${LCI_PKGLIST_X_LALAPPS} lalapps
then
    echo 'WARNING: no LALSuite release available' > /no-lalsuite
    cat /no-lalsuite
fi

# disable Python from downloading packages on-the-fly
# - everything should be specified as BuildRequires
# - requiring a virtual environment should disable pip install
python3 -m pip config set global.require-virtualenv True

# print info
dnf list installed --quiet
python3 -m pip list

# cleanup
dnf clean all
python3 -m pip cache purge
cd /tmp
rm -rf /tmp/*

EOF
