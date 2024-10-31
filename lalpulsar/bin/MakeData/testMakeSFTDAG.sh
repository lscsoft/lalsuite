## common variables
seg1_tstart=1257741529
seg1_tend=1257743329
seg2_tstart=1257743330
seg2_tend=1257745130
Tsft=1800
fmin=10
Band=1990

## a segment file
segs="./segs"
echo "${seg1_tstart} ${seg1_tend}" > $segs
echo "${seg2_tstart} ${seg2_tend}" >> $segs

## channels and expected filenames
pubopts="-O 4 -K DEV -R 1"
privopts="-O 0 -X private"
chan1="H1:GDS-CALIB_STRAIN_CLEAN"
chan2="H1:GDS-CALIB_STRAIN"
chan1sft1="H-1_H1_${Tsft}SFT_O4DEV+R1+CGDSCALIBSTRAINCLEAN+WHANN-${seg1_tstart}-${Tsft}.sft"
chan1sft2="H-1_H1_${Tsft}SFT_O4DEV+R1+CGDSCALIBSTRAINCLEAN+WHANN-${seg2_tstart}-${Tsft}.sft"
chan2sft1="H-1_H1_${Tsft}SFT_O4DEV+R1+CGDSCALIBSTRAIN+WHANN-${seg1_tstart}-${Tsft}.sft"
chan2sft2="H-1_H1_${Tsft}SFT_O4DEV+R1+CGDSCALIBSTRAIN+WHANN-${seg2_tstart}-${Tsft}.sft"
sft1priv="H-1_H1_${Tsft}SFT_private-${seg1_tstart}-1800.sft"
sft2priv="H-1_H1_${Tsft}SFT_private-${seg2_tstart}-1800.sft"

## public options
## run lalpulsar_MakeSFTDAG to create a fake output
cmdline="./lalpulsar_MakeSFTDAG ${pubopts} -f test.dag -G TEST -d H1_HOFT_C00 -k 7 -T ${Tsft} -p . -N ${chan1} -F ${fmin} -B ${Band} -w hann -P 0.5 -m 1 -A ligo.sim.o4.cw.explore.test -U albert.einstein -g segs -J /tmp/path/to"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
dagfile="test.dag"
datafindsub="datafind.sub"
sftsub="MakeSFTs.sub"
for file in $dagfile $datafindsub $sftsub; do
    if ! test -f $file; then
        echo "ERROR: could not find file '$file'"
        exit 1
    fi
done

## check dag content
testdagcontent=$(<$dagfile)
dagfilecontent="JOB datafind_1 datafind.sub
RETRY datafind_1 1
VARS datafind_1 gpsstarttime=\"1257741529\" gpsendtime=\"1257743329\" observatory=\"H\" inputdatatype=\"H1_HOFT_C00\" tagstring=\"TEST_1\"
JOB MakeSFTs_1 MakeSFTs.sub
RETRY MakeSFTs_1 1
VARS MakeSFTs_1 argList=\"${pubopts} -f 7 -t 1800 -p . -C H-1257741529-1257743329.cache -s 1257741529 -e 1257743329 -N ${chan1} -F 10 -B 1990 -w hann -P 0.5\" cachefile=\"cache/H-1257741529-1257743329.cache\" outputfiles=\"${chan1sft1}\" remapfiles=\"${chan1sft1}=${chan1sft1}\" tagstring=\"TEST_1\"
PARENT datafind_1 CHILD MakeSFTs_1
JOB datafind_2 datafind.sub
RETRY datafind_2 1
VARS datafind_2 gpsstarttime=\"1257743330\" gpsendtime=\"1257745130\" observatory=\"H\" inputdatatype=\"H1_HOFT_C00\" tagstring=\"TEST_2\"
JOB MakeSFTs_2 MakeSFTs.sub
RETRY MakeSFTs_2 1
VARS MakeSFTs_2 argList=\"${pubopts} -f 7 -t 1800 -p . -C H-1257743330-1257745130.cache -s 1257743330 -e 1257745130 -N ${chan1} -F 10 -B 1990 -w hann -P 0.5\" cachefile=\"cache/H-1257743330-1257745130.cache\" outputfiles=\"${chan1sft2}\" remapfiles=\"${chan1sft2}=${chan1sft2}\" tagstring=\"TEST_2\"
PARENT datafind_2 CHILD MakeSFTs_2"
if ! [[ $testdagcontent == $dagfilecontent ]]; then
   echo "ERROR: dagfile content did not match expected content"
   echo "test content:"
   echo $testdagcontent
   echo "Expected content:"
   echo $dagfilecontent
   exit 1
fi

## check submit files contents
testdatafindcontent=$(<$datafindsub)
datafindfilecontent="universe = vanilla
executable = /usr/bin/gw_data_find
arguments = --observatory \$(observatory) --url-type file --gps-start-time \$(gpsstarttime) --gps-end-time \$(gpsendtime) --lal-cache --gaps --type \$(inputdatatype)
getenv = *DATAFIND*
request_disk = 5MB
request_memory = 2000MB
accounting_group = ligo.sim.o4.cw.explore.test
accounting_group_user = albert.einstein
log = logs/datafind_test.dag.log
error = logs/datafind_\$(tagstring).err
output = cache/\$(observatory)-\$(gpsstarttime)-\$(gpsendtime).cache
should_transfer_files = yes
notification = never
queue 1"
if ! [[ $testdatafindcontent == $datafindfilecontent ]]; then
   echo "ERROR: datafind.sub content did not match expected content"
   echo "test content:"
   echo $testdatafindcontent
   echo "Expected content:"
   echo $datafindfilecontent
   exit 1
fi

testsftsubcontent=$(<$sftsub)
sftsubfilecontent="universe = vanilla
executable = /tmp/path/to/lalpulsar_MakeSFTs
arguments = \$(argList)
accounting_group = ligo.sim.o4.cw.explore.test
accounting_group_user = albert.einstein
log = logs/MakeSFTs_test.dag.log
error = logs/MakeSFTs_\$(tagstring).err
output = logs/MakeSFTs_\$(tagstring).out
notification = never
request_memory = 2048MB
request_disk = 1024MB
RequestCpus = 1
should_transfer_files = yes
transfer_input_files = \$(cachefile)
transfer_output_files = \$(outputfiles)
transfer_output_remaps = \"\$(remapfiles)\"
queue 1"
if ! [[ $testsftsubcontent == $sftsubfilecontent ]]; then
   echo "ERROR: MakeSFT.sub content did not match expected content"
   echo "test content:"
   echo $testsftsubcontent
   echo "Expected content:"
   echo $sftsubfilecontent
   exit 1
fi

## run lalpulsar_MakeSFTDAG to create a fake output for 2 channels
cmdline="./lalpulsar_MakeSFTDAG ${pubopts} -f test.dag -G TEST -d H1_HOFT_C00 -k 7 -T ${Tsft} -p dir1 dir2 -N ${chan1} ${chan2} -F ${fmin} -B ${Band} -w hann -P 0.5 -m 1 -A ligo.sim.o4.cw.explore.test -U albert.einstein -g segs -J /tmp/path/to"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
dagfile="test.dag"
datafindsub="datafind.sub"
sftsub="MakeSFTs.sub"
for file in $dagfile $datafindsub $sftsub; do
    if ! test -f $file; then
        echo "ERROR: could not find file '$file'"
        exit 1
    fi
done

## check dag content
testdagcontent=$(<$dagfile)
dagfilecontent="JOB datafind_1 datafind.sub
RETRY datafind_1 1
VARS datafind_1 gpsstarttime=\"1257741529\" gpsendtime=\"1257743329\" observatory=\"H\" inputdatatype=\"H1_HOFT_C00\" tagstring=\"TEST_1\"
JOB MakeSFTs_1 MakeSFTs.sub
RETRY MakeSFTs_1 1
VARS MakeSFTs_1 argList=\"${pubopts} -f 7 -t 1800 -p .,. -C H-1257741529-1257743329.cache -s 1257741529 -e 1257743329 -N ${chan1},${chan2} -F 10 -B 1990 -w hann -P 0.5\" cachefile=\"cache/H-1257741529-1257743329.cache\" outputfiles=\"${chan1sft1},${chan2sft1}\" remapfiles=\"${chan1sft1}=dir1/${chan1sft1};${chan2sft1}=dir2/${chan2sft1}\" tagstring=\"TEST_1\"
PARENT datafind_1 CHILD MakeSFTs_1
JOB datafind_2 datafind.sub
RETRY datafind_2 1
VARS datafind_2 gpsstarttime=\"1257743330\" gpsendtime=\"1257745130\" observatory=\"H\" inputdatatype=\"H1_HOFT_C00\" tagstring=\"TEST_2\"
JOB MakeSFTs_2 MakeSFTs.sub
RETRY MakeSFTs_2 1
VARS MakeSFTs_2 argList=\"${pubopts} -f 7 -t 1800 -p .,. -C H-1257743330-1257745130.cache -s 1257743330 -e 1257745130 -N ${chan1},${chan2} -F 10 -B 1990 -w hann -P 0.5\" cachefile=\"cache/H-1257743330-1257745130.cache\" outputfiles=\"${chan1sft2},${chan2sft2}\" remapfiles=\"${chan1sft2}=dir1/${chan1sft2};${chan2sft2}=dir2/${chan2sft2}\" tagstring=\"TEST_2\"
PARENT datafind_2 CHILD MakeSFTs_2"
if ! [[ $testdagcontent == $dagfilecontent ]]; then
   echo "ERROR: dagfile content did not match expected content"
   echo "test content:"
   echo $testdagcontent
   echo "Expected content:"
   echo $dagfilecontent
   exit 1
fi

## run lalpulsar_MakeSFTDAG to create a fake output for 1 channel private SFTs
cmdline="./lalpulsar_MakeSFTDAG ${privopts} -f test.dag -G TEST -d H1_HOFT_C00 -k 7 -T ${Tsft} -p . -N ${chan1} -F ${fmin} -B ${Band} -w hann -P 0.5 -m 1 -A ligo.sim.o4.cw.explore.test -U albert.einstein -g segs -J /tmp/path/to"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
dagfile="test.dag"
datafindsub="datafind.sub"
sftsub="MakeSFTs.sub"
for file in $dagfile $datafindsub $sftsub; do
    if ! test -f $file; then
        echo "ERROR: could not find file '$file'"
        exit 1
    fi
done

## check dag content
testdagcontent=$(<$dagfile)
dagfilecontent="JOB datafind_1 datafind.sub
RETRY datafind_1 1
VARS datafind_1 gpsstarttime=\"1257741529\" gpsendtime=\"1257743329\" observatory=\"H\" inputdatatype=\"H1_HOFT_C00\" tagstring=\"TEST_1\"
JOB MakeSFTs_1 MakeSFTs.sub
RETRY MakeSFTs_1 1
VARS MakeSFTs_1 argList=\"${privopts} -f 7 -t 1800 -p . -C H-1257741529-1257743329.cache -s 1257741529 -e 1257743329 -N ${chan1} -F 10 -B 1990 -w hann -P 0.5\" cachefile=\"cache/H-1257741529-1257743329.cache\" outputfiles=\"${sft1priv}\" remapfiles=\"${sft1priv}=${sft1priv}\" tagstring=\"TEST_1\"
PARENT datafind_1 CHILD MakeSFTs_1
JOB datafind_2 datafind.sub
RETRY datafind_2 1
VARS datafind_2 gpsstarttime=\"1257743330\" gpsendtime=\"1257745130\" observatory=\"H\" inputdatatype=\"H1_HOFT_C00\" tagstring=\"TEST_2\"
JOB MakeSFTs_2 MakeSFTs.sub
RETRY MakeSFTs_2 1
VARS MakeSFTs_2 argList=\"${privopts} -f 7 -t 1800 -p . -C H-1257743330-1257745130.cache -s 1257743330 -e 1257745130 -N ${chan1} -F 10 -B 1990 -w hann -P 0.5\" cachefile=\"cache/H-1257743330-1257745130.cache\" outputfiles=\"${sft2priv}\" remapfiles=\"${sft2priv}=${sft2priv}\" tagstring=\"TEST_2\"
PARENT datafind_2 CHILD MakeSFTs_2"
if ! [[ $testdagcontent == $dagfilecontent ]]; then
   echo "ERROR: dagfile content did not match expected content"
   echo "test content:"
   echo $testdagcontent
   echo "Expected content:"
   echo $dagfilecontent
   exit 1
fi

## run lalpulsar_MakeSFTDAG to create a fake output using a frame cache file
cmdline="./lalpulsar_MakeSFTDAG ${privopts} -f test.dag -G TEST -d H1_HOFT_C00 -k 7 -T ${Tsft} -p . -N ${chan1} -F ${fmin} -B ${Band} -w hann -P 0.5 -m 1 -A ligo.sim.o4.cw.explore.test -U albert.einstein -g segs -J /tmp/path/to -e /tmp/path/to.cache"
if ! eval "$cmdline"; then
    echo "ERROR: something failed when running '$cmdline'"
    exit 1
fi
dagfile="test.dag"
sftsub="MakeSFTs.sub"
for file in $dagfile $sftsub; do
    if ! test -f $file; then
        echo "ERROR: could not find file '$file'"
        exit 1
    fi
done

testdagcontent=$(<$dagfile)
dagfilecontent="JOB MakeSFTs_1 MakeSFTs.sub
RETRY MakeSFTs_1 1
VARS MakeSFTs_1 argList=\"${privopts} -f 7 -t 1800 -p . -C to.cache -s 1257741529 -e 1257743329 -N ${chan1} -F 10 -B 1990 -w hann -P 0.5\" cachefile=\"/tmp/path/to.cache\" outputfiles=\"${sft1priv}\" remapfiles=\"${sft1priv}=${sft1priv}\" tagstring=\"TEST_1\"
JOB MakeSFTs_2 MakeSFTs.sub
RETRY MakeSFTs_2 1
VARS MakeSFTs_2 argList=\"${privopts} -f 7 -t 1800 -p . -C to.cache -s 1257743330 -e 1257745130 -N ${chan1} -F 10 -B 1990 -w hann -P 0.5\" cachefile=\"/tmp/path/to.cache\" outputfiles=\"${sft2priv}\" remapfiles=\"${sft2priv}=${sft2priv}\" tagstring=\"TEST_2\""
if ! [[ $testdagcontent == $dagfilecontent ]]; then
   echo "ERROR: dagfile content did not match expected content"
   echo "test content:"
   echo $testdagcontent
   echo "Expected content:"
   echo $dagfilecontent
   exit 1
fi

rm MakeSFTs.sub datafind.sub test.dag segs
rmdir logs cache
