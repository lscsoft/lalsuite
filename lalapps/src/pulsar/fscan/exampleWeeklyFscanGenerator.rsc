# exampleFscanGenerator.rsc

# Set tcl list with list of channels, frame types, etc....
# Each item in the list contains: channel_name channel_type IFO input_type output_sft_dir comparison_chan comparison_snr comparison_delta_dir knee_freq sft_time_baseline start_freq band sub_band extra_time_for_data_find freq_res alternative_dir_name
# If comparison_delta_dir is 0 compare with current fscans, if -1 compare with last previous set of fscans, etc...
set ::masterList {\
{H2:LSC-STRAIN ADC_REAL8 H2 H2_RDS_C03_L2 default H2:LSC-STRAIN 4 -1 30 1800 50 951 50 64 0.1 default}\
{H0:PEM-BSC9_MAGX ADC_INT2 H2 RDS_R_L1 default H2:LSC-STRAIN 4 0 30 1800 50 951 50 64 0.1 default}\
}

# If fixedComparison is 1 then used then ignore the comparison delta dir in the masterList but compare using fixd values: 
set fixedComparison 1;
set fixedComparisonChanDir ../../fscans_2009_03_30_13_02_17_PDT_Mon/H2_LSC-STRAIN;
set fixedComparisonString "H2_923700000_923714400";
set fixedComparisonSNR  4;

set fscanDriverPath "/home/gmendell/bin/fscanDriver.py"; # complete path to fscanDriver.py
set matlabPath "/ldcg/matlab_r2008a";      # Path to matlab installation to use with -m option to fscanDriver.py, e.g., /ldcg/matlab_r2008a

set parentOutputDirectory "weekly";

set moveSFTsFromDirList "/usr1/ugmendell/oldsearchcode/src/lalapps/src/pulsar/fscan/testDaily/daily";
set moveSFTsFromSuffix "sfts/tmp";

set startTime 923628000;
#set endTime "now";
# Time lag before now to end fscans, if endTime set to now: 
set timeLag 7200;
set endTime "935448588";
set startTimeFile "lastTimeUsedByFscanGenerator.txt";

set useEndTimeForDirName 1;    # 0 == false, 1 == true; Use start time to name output directory.

# Control whether to find segs using ligolw_segment_query:
set useLSCsegFind 0 
set typeLSCsegFind "H2:DMT-SCIENCE:1"
#set ::segFindCmdAndPath /opt/lscsoft/glue/bin/ligolw_segment_query
set ::segFindCmdAndPath /bin/echo 
set ::grepCommandAndPath /bin/grep
#set ::ligolwPrintCommandAndPath /opt/lscsoft/glue/bin/ligolw_print
set ::ligolwPrintCommandAndPath /bin/echo 

