FROM igwn/base:buster

LABEL name="LALSuite Nightly - Debian Buster" \
      maintainer="Adam Mercer <adam.mercer@ligo.org>" \
      support="Not Supported"

# add debian packages to container
COPY debs /srv/local-apt=repository

# install debs & cleanup
RUN apt-get update && \
      apt-get -y install local-apt-repository && \
      /usr/lib/local-apt-repository/rebuild && \
      apt-get update && \
      apt-get install lscsoft-lalsuite && \
      rm -rf /var/lib/apts/lists/*
