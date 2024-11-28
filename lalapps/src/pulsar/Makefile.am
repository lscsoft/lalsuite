MOSTLYCLEANFILES =
dist_bin_SCRIPTS =
include $(top_srcdir)/gnuscripts/lalsuite_old_exec_name_stubs.am

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

CLEANFILES = \
	CreateEphemeris/lalapps_create_time_correction_ephemeris \
	CrossCorr/lalapps_combine_crosscorr_toplists \
	CrossCorr/lalapps_pulsar_crosscorr \
	CrossCorr/lalapps_pulsar_crosscorr_v2 \
	CrossCorr/lalapps_run_pulsar_crosscorr_v2 \
	FITSTools/lalapps_fits_header_getval \
	FITSTools/lalapps_fits_header_list \
	FITSTools/lalapps_fits_overview \
	FITSTools/lalapps_fits_table_list \
	Fscan/lalapps_spec_avg \
	Fscan/lalapps_spec_avg_long \
	Fscan/lalapps_spec_coherence \
	Fstatistic/lalapps_compareFstats \
	Fstatistic/lalapps_ComputeFstatBenchmark \
	Fstatistic/lalapps_ComputeFstatistic_v2 \
	Fstatistic/lalapps_ComputeFstatLatticeCount \
	Fstatistic/lalapps_ComputeFstatMCUpperLimit \
	Fstatistic/lalapps_PredictFstat \
	Fstatistic/lalapps_synthesizeBstatMC \
	Fstatistic/lalapps_synthesizeLVStats \
	Fstatistic/lalapps_synthesizeTransientStats \
	Fstatistic/SemiAnalyticF \
	GCT/lalapps_HierarchSearchGCT \
	HeterodyneSearch/lalapps_create_pulsar_signal_frame \
	HeterodyneSearch/lalapps_heterodyne_pulsar \
	HeterodyneSearch/lalapps_knope \
	HeterodyneSearch/lalapps_knope_automation_script \
	HeterodyneSearch/lalapps_knope_collate_results \
	HeterodyneSearch/lalapps_knope_result_page \
	HeterodyneSearch/lalapps_pulsar_frequency_evolution \
	HeterodyneSearch/lalapps_pulsar_parameter_estimation_nested \
	HeterodyneSearch/lalapps_SplInter \
	HeterodyneSearch/lalapps_ssbtodetector \
	HoughFstat/lalapps_HierarchicalSearch \
	Hough/lalapps_DriveHoughMulti \
	Hough/lalapps_HoughValidate \
	Hough/lalapps_HoughValidateAM \
	Hough/lalapps_MCInjectHoughMulti \
	Hough/lalapps_MCInjectHoughMultiChi2Test \
	Hough/lalapps_MultiWeights \
	Hough/lalapps_ValidateChi2Test \
	Hough/lalapps_ValidateHoughMulti \
	MakeData/compareTS \
	MakeData/lalapps_Makefakedata_v4 \
	MakeData/lalapps_Makefakedata_v5 \
	MakeData/lalapps_MakeSFTDAG \
	MakeData/lalapps_MakeSFTs \
	SFTTools/lalapps_compareSFTs \
	SFTTools/lalapps_ComputePSD \
	SFTTools/lalapps_dumpSFT \
	SFTTools/lalapps_SFTclean \
	SFTTools/lalapps_SFTvalidate \
	SFTTools/lalapps_splitSFTs \
	SFTTools/lalapps_WriteSFTsfromSFDBs \
	SFTTools/SFTwrite \
	SidebandSearch/lalapps_CombSearch \
	Tools/lalapps_ComputeAntennaPattern \
	Tools/lalapps_FstatMetric_v2 \
	Tools/lalapps_PrintDetectorState \
	TwoSpect/lalapps_TwoSpect \
	TwoSpect/lalapps_TwoSpectTemplateBank \
	Weave/lalapps_Weave \
	Weave/lalapps_WeaveCompare \
	Weave/lalapps_WeaveConcat \
	Weave/lalapps_WeaveSetup \
	$(END_OF_LIST)

clean-local:
	-rm -rf */*/.deps */*/.libs */.deps */.libs */*.testdir */__pycache__ */.pytest_cache
	-rm -f */.*.log */*.log */*.o */*.trs */*.xml */Makefile */Makefile.in
