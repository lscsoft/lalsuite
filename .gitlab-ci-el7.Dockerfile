FROM igwn/base:el7-testing

LABEL name="LALSuite Nightly - EL7" \
      maintainer="Adam Mercer <adam.mercer@ligo.org>" \
      support="Reference Platform"

# add RPMs to container
COPY rpms /rpms

# install RPMs & cleanup
RUN yum makecache && \
      yum -y localinstall /rpms/*.rpm && \
      rm -rf /rpms && yum clean all
