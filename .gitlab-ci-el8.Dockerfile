FROM igwn/base:el8-testing

LABEL name="LALSuite Nightly - EL8" \
      maintainer="Adam Mercer <adam.mercer@ligo.org>" \
      support="Reference Platform"

# add RPMs to container
COPY rpms /rpms

# install latest release
RUN dnf -y -q upgrade && \
      dnf -y -q install createrepo && \
      createrepo --quiet /rpms && \
      dnf clean all

# generate repo metadata
RUN echo "[local-builds]" > /etc/yum.repos.d/local-builds.repo && \
    echo "name=Local builds" >> /etc/yum.repos.d/local-builds.repo && \
    echo "baseurl=file:///rpms" >> /etc/yum.repos.d/local-builds.repo && \
    echo "enabled=1" >> /etc/yum.repos.d/local-builds.repo && \
    echo "gpgcheck=0" >> /etc/yum.repos.d/local-builds.repo

# install packages from local repo
RUN dnf -y install lal* python*-lal --exclude lalsuite* && \
      dnf clean all
