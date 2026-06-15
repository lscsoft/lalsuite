# syntax=docker/dockerfile:1

ARG EL_VERSION

FROM igwn/base:${EL_VERSION}

ARG EL_VERSION

LABEL name="LALSuite Release Image - Enterprise Linux ${EL_VERSION}"
LABEL maintainer="LALSuite Maintainers <lal-discuss@ligo.org>"
LABEL support="Reference Platform"

SHELL ["/bin/bash", "-c"]

# copy packages
COPY ./packages/ /packages/

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

# install required packages
dnf -y -q install \
    findutils \
    ;

# install upstream RPMs
upstream_rpms=$(find /packages -name '*.rpm')
echo "===== upstream RPMs"
printf "%s\n" ${upstream_rpms}
echo "====="
dnf install -y ${upstream_rpms}

# print info
dnf list installed --quiet

# cleanup
dnf clean all
rm -rf /packages

EOF
