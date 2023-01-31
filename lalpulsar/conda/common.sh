# common build options for LALPulsar

# macros
_make="make -j ${CPU_COUNT} V=1 VERBOSE=1"

# -- compile customisations

# replace package name in debug-prefix-map with source name
export CFLAGS=$(
   echo ${CFLAGS:-} |
   sed -E 's|\/usr\/local\/src\/conda\/'${PKG_NAME}'|/usr/local/src/conda/lalpulsar|g'
)

# link options
export GSL_LIBS="-L${PREFIX}/lib -lgsl"

# -- configure arguments

CONFIGURE_ARGS="
  --disable-doxygen
  --disable-static
  --disable-swig-octave
  --enable-cfitsio
  --enable-openmp
  --prefix=${PREFIX}
"

# customisation for LALSuite development CI
if [[ "${GITLAB_CI}" == "true" ]] && [[ "x${CI_COMMIT_TAG}" == x ]]; then
  # declare nightly builds
  if [ ! -z "${ENABLE_NIGHTLY}" ]; then
    CONFIGURE_ARGS="${CONFIGURE_ARGS} --enable-nightly"
  fi
# production builds ignore GCC warnings
else
  CONFIGURE_ARGS="${CONFIGURE_ARGS} --disable-gcc-flags"
fi

# disable help2man when cross-compiling
if [[ "${CONDA_BUILD_CROSS_COMPILATION:-}" = "1" && "${CROSSCOMPILING_EMULATOR}" = "" ]]; then
  CONFIGURE_ARGS="${CONFIGURE_ARGS} --disable-help2man"
fi
