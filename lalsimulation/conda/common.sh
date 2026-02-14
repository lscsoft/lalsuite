# common build options for LALSimulation

# macros
_make="make -j ${CPU_COUNT} V=1 VERBOSE=1"

# -- compile customisations

# replace package name in debug-prefix-map with source name
export CFLAGS=$(
   echo ${CFLAGS:-} |
   sed -E 's|\/usr\/local\/src\/conda\/'${PKG_NAME}'|/usr/local/src/conda/lalsimulation|g'
)

# link options
export GSL_LIBS="-L${PREFIX}/lib -lgsl"

# -- configure arguments

CONFIGURE_ARGS="
  --disable-doxygen
  --disable-static
  --disable-swig-octave
  --enable-openmp
  --prefix=${PREFIX}
"

# customisation for LALSuite development CI
if [[ "${GITLAB_CI}" == "true" ]]; then
  # declare nightly builds
  CONFIGURE_ARGS="${CONFIGURE_ARGS} ${ENABLE_NIGHTLY}"
else
  # production builds ignore GCC warnings
  CONFIGURE_ARGS="${CONFIGURE_ARGS} --disable-gcc-flags"
fi

# disable help2man when cross-compiling
if [[ "${CONDA_BUILD_CROSS_COMPILATION:-}" = "1" && "${CROSSCOMPILING_EMULATOR}" = "" ]]; then
  CONFIGURE_ARGS="${CONFIGURE_ARGS} --disable-help2man"
fi
