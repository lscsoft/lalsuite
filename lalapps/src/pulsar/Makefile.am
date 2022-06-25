MOSTLYCLEANFILES =
dist_bin_SCRIPTS =
include $(top_srcdir)/gnuscripts/lalsuite_old_exec_name_stubs.am

SUBDIRS = \
	Tools \
	SFTTools \
	FITSTools \
	MakeData \
	Fscan \
	CreateEphemeris \
	Fstatistic \
	HeterodyneSearch \
	Hough \
	HoughFstat \
	GCT \
	CrossCorr \
	SidebandSearch \
	TwoSpect \
	Weave \
	$(END_OF_LINE)

# Because many tests in lalapps/src/pulsar/ call executables from other
# subdirectories in lalapps/src/pulsar/, it is safest to make sure the
# whole of lalapps/src/pulsar/ is built first
check-recursive: all-recursive

# All programs under lalapps/src/pulsar/ [1] must be executed by at least one
# test script run by "make check". This applies to C executables, Python
# scripts, etc. The following rule enforces this requirement.
#
# [1] Certain directories are currently exempt from this requirement; in these
# directories, Makefile.am includes lalapps.am instead of lalapps_pulsar.am
.PHONY: clean-execs-tested-logs check-execs-are-tested

check-recursive: clean-execs-tested-logs

check-local: check-execs-are-tested

clean-execs-tested-logs:
	@tested_logfiles=`find $(builddir) -name '.tested_*.log'`; \
	for file in $${tested_logfiles}; do \
		echo '# empty' > "$${file}"; \
	done

check-execs-are-tested:
	@tested_logfiles=`find $(builddir) -name '.tested_*.log'`; \
	not_tested=`$(AWK) ' \
	BEGIN { anyprog = 0 } \
	$$1 == "programs" { anyprog = 1; for (i = 3; i <= NF; ++i) { prog[$$i] = 1; progpath[$$i] = $$2 "/" $$i } } \
	$$1 == "tested" { for (i = 2; i <= NF; ++i) { tested[$$i] = 1; } } \
	END { if (anyprog) { for (e in prog) { if (!tested[e]) { print progpath[e]; } } } else { print "NONE"; } } \
	' $${tested_logfiles} /dev/null | $(SED) 's|^$(subdir)/||' | sort`; \
	if test "X$${not_tested}" = XNONE; then \
		exit 0; \
	fi; \
	if test "X$${not_tested}" != X; then \
		printf "\nERROR: the following executables under lalapps_pulsar.am test harness are not tested:\n" >&2; \
		printf "   %s\n" $${not_tested} >&2; \
		printf "\n" >&2; \
		exit 1; \
	fi; \
	printf "OK: all C executables under lalapps/src/pulsar/ have been tested\n"

old_exec_name_stubs = \
	lalapps_ComputeAntennaPattern \
	lalapps_ComputeFstatBenchmark \
	lalapps_ComputeFstatLatticeCount \
	lalapps_ComputeFstatMCUpperLimit \
	lalapps_ComputeFstatistic_v2 \
	lalapps_ComputePSD \
	lalapps_CopySFTs \
	lalapps_DriveHoughMulti \
	lalapps_FstatMetric_v2 \
	lalapps_HierarchSearchGCT \
	lalapps_HierarchicalSearch \
	lalapps_MakeSFTDAG \
	lalapps_MakeSFTs \
	lalapps_Makefakedata_v4 \
	lalapps_Makefakedata_v5 \
	lalapps_PredictFstat \
	lalapps_PrintDetectorState \
	lalapps_SFTclean \
	lalapps_SFTvalidate \
	lalapps_Weave \
	lalapps_WeaveCompare \
	lalapps_WeaveConcat \
	lalapps_WeaveSetup \
	lalapps_WriteSFTsfromSFDBs \
	lalapps_combine_crosscorr_toplists \
	lalapps_compareFstats \
	lalapps_compareSFTs \
	lalapps_create_solar_system_ephemeris \
	lalapps_create_solar_system_ephemeris_python \
	lalapps_create_time_correction_ephemeris \
	lalapps_dumpSFT \
	lalapps_fits_header_getval \
	lalapps_fits_header_list \
	lalapps_fits_overview \
	lalapps_fits_table_list \
	lalapps_heterodyne_pulsar \
	lalapps_knope \
	lalapps_knope_automation_script \
	lalapps_knope_collate_results \
	lalapps_knope_result_page \
	lalapps_pulsar_crosscorr_v2 \
	lalapps_pulsar_frequency_evolution \
	lalapps_pulsar_parameter_estimation_nested \
	lalapps_run_pulsar_crosscorr_v2 \
	lalapps_spec_avg \
	lalapps_spec_avg_long \
	lalapps_spec_coherence \
	lalapps_splitSFTs \
	lalapps_ssbtodetector \
	lalapps_synthesizeBstatMC \
	lalapps_synthesizeLVStats \
	lalapps_synthesizeTransientStats \
	$(END_OF_LINE)

new_exec_name_regex = ' \
	s|^lalapps_pulsar_crosscorr_v2$$|lalpulsar_crosscorr_v2|; \
	s|^lalapps_pulsar_frequency_evolution$$|lalpulsar_frequency_evolution|; \
	s|^lalapps_heterodyne_pulsar$$|lalpulsar_heterodyne|; \
	s|^lalapps_pulsar_parameter_estimation_nested$$|lalpulsar_parameter_estimation_nested|; \
	s|^lalapps_run_pulsar_crosscorr_v2$$|lalpulsar_run_crosscorr_v2|; \
	s|^lalapps_|lalpulsar_|'

EXTRA_DIST = old_exec_name_stubs.dox
