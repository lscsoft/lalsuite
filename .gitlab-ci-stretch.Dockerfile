FROM igwn/base:stretch-proposed

LABEL name="LALSuite Nightly - Debian Stretch" \
      maintainer="Adam Mercer <adam.mercer@ligo.org>" \
      support="Not Supported"

# add debian packages to container
COPY debs /debs

# install debs & cleanup
RUN apt-get update && \
      dpkg -i /debs/*.deb || apt-get --assume-yes -f install && \
      rm -rf /debs && rm -rf /var/lib/apts/lists/*
