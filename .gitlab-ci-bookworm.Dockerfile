FROM igwn/base:bookworm

LABEL name="LALSuite Nightly - Debian Bullseye" \
      maintainer="Adam Mercer <adam.mercer@ligo.org>" \
      support="Not Supported"

# add debian packages to container
COPY debs /srv/local-apt-repository

# install debs & cleanup
RUN apt-get update && \
      apt-get -y upgrade && \
      apt-get -y install lalsuite lalsuite-dev && \
      dpkg -P lalsuite lalsuite-dev && \
      apt-get -y install local-apt-repository && \
      /usr/lib/local-apt-repository/rebuild && \
      apt-get update && \
      apt-get -y upgrade && \
      rm -rf /var/lib/apts/lists/*
