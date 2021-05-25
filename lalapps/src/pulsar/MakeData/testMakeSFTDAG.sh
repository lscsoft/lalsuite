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

## run lalapps_MakeSFTDAG to create a fake output
cmdline="lalapps_MakeSFTDAG -f test.dag -G TEST -d H1_HOFT_C00 -k 7 -T ${Tsft} -p . -N H1:GDS-CALIB_STRAIN_CLEAN -F ${fmin} -B ${Band} -D 3 -X TEST -w 3 -P 0.5 -m 1 -A ligo.sim.o4.cw.explore.test -U albert.einstein -H -g segs"
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

testdagcontent=$(<$dagfile)
dagfilecontent="JOB LSCdataFind_1 datafind.sub
RETRY LSCdataFind_1 10
VARS LSCdataFind_1 gpsstarttime=\"1257741529\" gpsendtime=\"1257743329\" observatory=\"H\" inputdatatype=\"H1_HOFT_C00\" tagstring=\"TEST_1\"
JOB MakeSFTs_1 MakeSFTs.sub
RETRY MakeSFTs_1 5
VARS MakeSFTs_1 argList=\"-f 7 -t 1800 -p . -C cache/H-1257741529-1257743329.cache -s 1257741529 -e 1257743329 -N H1:GDS-CALIB_STRAIN_CLEAN -v 2 -F 10 -B 1990 -D 3 -X TEST -w 3 -P 0.5 -H\" tagstring=\"TEST_1\"
PARENT LSCdataFind_1 CHILD MakeSFTs_1
JOB LSCdataFind_2 datafind.sub
RETRY LSCdataFind_2 10
VARS LSCdataFind_2 gpsstarttime=\"1257743330\" gpsendtime=\"1257745130\" observatory=\"H\" inputdatatype=\"H1_HOFT_C00\" tagstring=\"TEST_2\"
JOB MakeSFTs_2 MakeSFTs.sub
RETRY MakeSFTs_2 5
VARS MakeSFTs_2 argList=\"-f 7 -t 1800 -p . -C cache/H-1257743330-1257745130.cache -s 1257743330 -e 1257745130 -N H1:GDS-CALIB_STRAIN_CLEAN -v 2 -F 10 -B 1990 -D 3 -X TEST -w 3 -P 0.5 -H\" tagstring=\"TEST_2\"
PARENT LSCdataFind_2 CHILD MakeSFTs_2"
if ! [[ $testdagcontent == $dagfilecontent ]]; then
   echo "ERROR: dagfile content did not match expected content"
   echo "test content:"
   echo $testdagcontent
   echo "Expected content:"
   echo $dagfilecontent
   exit 1
fi

testdatafindcontent=$(<$datafindsub)
datafindfilecontent="universe = vanilla
executable = gw_data_find
arguments = -r \$ENV(LIGO_DATAFIND_SERVER) --observatory \$(observatory) --url-type file --gps-start-time \$(gpsstarttime) --gps-end-time \$(gpsendtime) --lal-cache --type \$(inputdatatype) 
getenv = True
accounting_group = ligo.sim.o4.cw.explore.test
accounting_group_user = albert.einstein
log = logs/datafind_test.dag.log
error = logs/datafind_\$(tagstring).err
output = cache/\$(observatory)-\$(gpsstarttime)-\$(gpsendtime).cache
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
executable = lalapps_MakeSFTs
arguments = \$(argList)
getenv = True
accounting_group = ligo.sim.o4.cw.explore.test
accounting_group_user = albert.einstein
log = logs/MakeSFTs_test.dag.log
error = logs/MakeSFTs_\$(tagstring).err
output = logs/MakeSFTs_\$(tagstring).out
notification = never
RequestMemory = 2048
RequestCpus = 1
queue 1"
if ! [[ $testsftsubcontent == $sftsubfilecontent ]]; then
   echo "ERROR: MakeSFT.sub content did not match expected content"
   echo "test content:"
   echo $testsftsubcontent
   echo "Expected content:"
   echo $sftsubfilecontent
   exit 1
fi

rm MakeSFTs.sub datafind.sub segs test.dag
rmdir cache logs
