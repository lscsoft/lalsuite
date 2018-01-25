FROM ligo/lalsuite-runtime:jessie
COPY /opt/lalsuite /opt/lalsuite
ENV LD_LIBRARY_PATH="/opt/lalsuite/lib" \
    OCTAVE_PATH="/opt/lalsuite/lib/x86_64-linux-gnu/octave/3.8.2/site/oct/x86_64-pc-linux-gnu" \
    PATH="/opt/lalsuite/bin:${PATH}" \
    PKG_CONFIG_PATH="/opt/lalsuite/lib/pkgconfig" \
    PYTHONPATH="/opt/lalsuite/lib/python2.7/site-packages"

RUN find /opt/lalsuite/lib/pkgconfig -type f -name 'lal*.pc' -exec sed -i 's|/builds/lscsoft/lalsuite||g' {} \;

COPY /docker/startup /usr/local/bin/startup

ENTRYPOINT [ "/usr/local/bin/startup" ]
CMD ["/bin/bash", "-i" ]
