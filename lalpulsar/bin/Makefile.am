BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
TESTS =
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

TESTS += check_execs_are_tested
EXTRA_DIST += check_execs_are_tested

TESTS_ENVIRONMENT = \
	$(TESTS_ENV_CONFIG_VARS) \
	export builddir; builddir='$(builddir)';
	$(END_OF_LIST)

# -- C programs -------------

bin_PROGRAMS = \
	lalpulsar_version \
	$(END_OF_LIST)

lalpulsar_version_SOURCES = version.c

TESTS += \
	lalpulsar_version \
	$(END_OF_LIST)

# -- Python programs --------

if HAVE_PYTHON

pybin_scripts = \
	$(END_OF_LIST)

TESTS += \
	$(END_OF_LIST)

endif

# -- Manuals ----------------

man1_MANS = $(help2man_MANS)
