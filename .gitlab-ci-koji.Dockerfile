FROM igwn/base:el7-testing

# install release RPMs
RUN yum -y install lalsuite && \
    yum clean all

# setup lalsuite-nightly repository
RUN echo "[lalsuite-nightly]" > /etc/yum.repos.d/lalsuite-nightly.repo && \
    echo "name = lalsuite-nightly" >> /etc/yum.repos.d/lalsuite-nightly.repo && \
    echo "baseurl = https://koji.ligo-la.caltech.edu/kojifiles/repos-dist/epel7-lalsuite/latest/\$basearch/" >> /etc/yum.repos.d/lalsuite-nightly.repo && \
    echo "enabled = 1" >> /etc/yum.repos.d/lalsuite-nightly.repo && \
    echo "gpgcheck = 0" >> /etc/yum.repos.d/lalsuite-nightly.repo

# update LALSuite from nightly repository
RUN yum -y remove lalsuite && \
    yum -y update && \
    yum clean all