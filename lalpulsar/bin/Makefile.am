BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
man1_MANS =
include $(top_srcdir)/gnuscripts/lalsuite_test.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

# -- Subdirectories -------------

SUBDIRS = \
	Tools \
	SFTTools \
	FITSTools \
	MakeData \
	Fscan \
	CreateEphemeris \
	HWInjection \
	Fstatistic \
	HeterodyneSearch \
	Hough \
	HoughFstat \
	GCT \
	CrossCorr \
	SidebandSearch \
	TwoSpect \
	Weave \
	PiecewiseSearch \
	longtests \
	$(END_OF_LINE)

# Because many tests in lalpulsar/bin/ call executables from other
# subdirectories in lalpulsar/bin/, it is safest to make sure the whole of
# lalpulsar/bin/ is built first
check-recursive: all-recursive

# All C executables under lalpulsar/bin/ must be executed by at least one test
# script run by "make check". The following rules enforce this requirement.
.PHONY: clean-execs-tested-logs
check-recursive: clean-execs-tested-logs
clean-execs-tested-logs:
	@tested_logfiles=`find $(builddir) -name '.tested_*.log'`; \
	for file in $${tested_logfiles}; do \
		echo '# empty' > "$${file}"; \
	done

test_scripts += check_execs_are_tested.sh

# -- C programs -------------

bin_PROGRAMS = \
	lalpulsar_version \
	$(END_OF_LIST)

lalpulsar_version_SOURCES = version.c

# Add compiled test programs to this variable
test_programs += lalpulsar_version

# Add shell, Python, etc. test scripts to this variable
test_scripts +=

# Add any helper programs required by tests to this variable
test_helpers +=

# -- Python programs --------

if HAVE_PYTHON

pybin_scripts = \
	$(END_OF_LIST)

# Add compiled test programs to this variable
test_programs +=

# Add shell, Python, etc. test scripts to this variable
test_scripts +=

# Add any helper programs required by tests to this variable
test_helpers +=

endif
