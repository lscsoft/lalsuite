FROM igwn/base:el7-testing

LABEL name="LALSuite Nightly - EL7" \
      maintainer="Adam Mercer <adam.mercer@ligo.org>" \
      support="Reference Platform"

# add RPMs to container
COPY rpms /rpms

# install latest release
RUN yum -y -q upgrade && \
      yum -y -q install createrepo && \
      createrepo --quiet /rpms && \
      yum clean all

# generate repo metadata
RUN echo "[local-builds]" > /etc/yum.repos.d/local-builds.repo && \
    echo "name=Local builds" >> /etc/yum.repos.d/local-builds.repo && \
    echo "baseurl=file:///rpms" >> /etc/yum.repos.d/local-builds.repo && \
    echo "enabled=1" >> /etc/yum.repos.d/local-builds.repo && \
    echo "gpgcheck=0" >> /etc/yum.repos.d/local-builds.repo

# install packages from local repo
RUN yum -y install lal* python*-lal --exclude lalsuite* && \
      yum clean all
