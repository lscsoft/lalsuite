# common build options for LAL
#
# to simplify repeatedly specifying things to build the various components
#

# macros
_make="make -j ${CPU_COUNT} V=1 VERBOSE=1"

# -- compile customisations

# replace package name in debug-prefix-map with source name
export CFLAGS=$(
   echo ${CFLAGS:-} |
   sed -E 's|\/usr\/local\/src\/conda\/'${PKG_NAME}'|/usr/local/src/conda/lal|g'
)

# only link libraries we actually use
export GSL_LIBS="-L${PREFIX}/lib -lgsl"
export HDF5_LIBS="-L${PREFIX}/lib -lhdf5 -lhdf5_hl"

# -- configure arguments

CONFIGURE_ARGS="
  --prefix=${PREFIX}
  --disable-doxygen
  --disable-static
  --disable-swig-octave
  --with-hdf5=yes
"

# customisation for LALSuite development CI
if [[ "${GITLAB_CI}" == "true" ]]; then
  # use configuration flags set by CI
  CONFIGURE_ARGS="${CONFIGURE_ARGS} ${LCI_CONFIGURE_ARGS}"
else
  # production builds ignore GCC warnings
  CONFIGURE_ARGS="${CONFIGURE_ARGS} --disable-gcc-flags"
fi

# disable help2man when cross-compiling
if [[ "${CONDA_BUILD_CROSS_COMPILATION:-}" = "1" && "${CROSSCOMPILING_EMULATOR}" = "" ]]; then
  CONFIGURE_ARGS="${CONFIGURE_ARGS} --disable-help2man"
fi

# select FFT implementation
if [[ "${fft_impl}" == "mkl" ]]; then
    CONFIGURE_ARGS="${CONFIGURE_ARGS} --enable-intelfft"
fi
