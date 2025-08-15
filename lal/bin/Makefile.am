BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
man1_MANS =
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_test.am

# -- C programs -------------

bin_PROGRAMS = \
	lal_cache \
	lal_simd_detect \
	lal_tconvert \
	lal_version \
	$(END_OF_LIST)

lal_cache_SOURCES = cache.c
lal_simd_detect_SOURCES = simd_detect.c
lal_tconvert_SOURCES = tconvert.c getdate.c getdate.h
lal_version_SOURCES = version.c

if FFTW3

bin_PROGRAMS += \
	lal_fftwf_wisdom \
	lal_fftw_wisdom \
	$(END_OF_LIST)

lal_fftwf_wisdom_SOURCES = fftwf_wisdom.c
lal_fftw_wisdom_SOURCES = fftw_wisdom.c

endif

# Add compiled test programs to this variable
test_programs += lal_simd_detect
test_programs += lal_version

# Add shell, Python, etc. test scripts to this variable
test_scripts += test_cache.sh
test_scripts += test_tconvert.sh

# Add any helper programs required by tests to this variable
test_helpers +=

MOSTLYCLEANFILES += \
	Z-ilwd_test_frame-600000000-1.gwf \
	$(END_OF_LIST)

EXTRA_DIST += \
	getdate.y \
	$(END_OF_LIST)

# -- Python programs --------

if HAVE_PYTHON

pybin_scripts = \
	lal_path2cache \
	lal_searchsum2cache \
	$(END_OF_LIST)

# Add shell, Python, etc. test scripts to this variable
test_scripts += test_path2cache.sh

installcheck-local:
	$(bindir)/lal_path2cache --help

endif
