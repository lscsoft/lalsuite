FROM ligo/base:stretch

LABEL name="LALSuite Nightly - Debian Stretch" \
      maintainer="Adam Mercer <adam.mercer@ligo.org>" \
      support="Reference Platform" \
      date="20180506"

# add debian packages to container
COPY debs /debs

# install debs & cleanup
RUN apt-get update && \
      dpkg -i /debs/*.deb || apt-get --assume-yes -f install && \
      rm -rf /debs && rm -rf /var/lib/apts/lists/*
