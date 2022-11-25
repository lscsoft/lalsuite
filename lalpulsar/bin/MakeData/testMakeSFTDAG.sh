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

for SFT_name_opts in "-O 4 -K DEV -R 1" "-O 0 -X private"; do

## run lalpulsar_MakeSFTDAG to create a fake output
cmdline="lalpulsar_MakeSFTDAG ${SFT_name_opts} -f test.dag -G TEST -d H1_HOFT_C00 -k 7 -T ${Tsft} -p . -N H1:GDS-CALIB_STRAIN_CLEAN -F ${fmin} -B ${Band} -w 3 -P 0.5 -m 1 -A ligo.sim.o4.cw.explore.test -U albert.einstein -g segs -J /tmp/path/to"
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
dagfilecontent="JOB datafind_1 datafind.sub
RETRY datafind_1 10
VARS datafind_1 gpsstarttime=\"1257741529\" gpsendtime=\"1257743329\" observatory=\"H\" inputdatatype=\"H1_HOFT_C00\" tagstring=\"TEST_1\"
JOB MakeSFTs_1 MakeSFTs.sub
RETRY MakeSFTs_1 5
VARS MakeSFTs_1 argList=\"${SFT_name_opts} -f 7 -t 1800 -p . -C cache/H-1257741529-1257743329.cache -s 1257741529 -e 1257743329 -N H1:GDS-CALIB_STRAIN_CLEAN -F 10 -B 1990 -w 3 -P 0.5\" tagstring=\"TEST_1\"
PARENT datafind_1 CHILD MakeSFTs_1
JOB datafind_2 datafind.sub
RETRY datafind_2 10
VARS datafind_2 gpsstarttime=\"1257743330\" gpsendtime=\"1257745130\" observatory=\"H\" inputdatatype=\"H1_HOFT_C00\" tagstring=\"TEST_2\"
JOB MakeSFTs_2 MakeSFTs.sub
RETRY MakeSFTs_2 5
VARS MakeSFTs_2 argList=\"${SFT_name_opts} -f 7 -t 1800 -p . -C cache/H-1257743330-1257745130.cache -s 1257743330 -e 1257745130 -N H1:GDS-CALIB_STRAIN_CLEAN -F 10 -B 1990 -w 3 -P 0.5\" tagstring=\"TEST_2\"
PARENT datafind_2 CHILD MakeSFTs_2"
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
executable = /usr/bin/gw_data_find
arguments = -r \$ENV(LIGO_DATAFIND_SERVER) --observatory \$(observatory) --url-type file --gps-start-time \$(gpsstarttime) --gps-end-time \$(gpsendtime) --lal-cache --type \$(inputdatatype) 
getenv = True
request_disk = 5MB
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
executable = /tmp/path/to/lalpulsar_MakeSFTs
arguments = \$(argList)
getenv = True
accounting_group = ligo.sim.o4.cw.explore.test
accounting_group_user = albert.einstein
log = logs/MakeSFTs_test.dag.log
error = logs/MakeSFTs_\$(tagstring).err
output = logs/MakeSFTs_\$(tagstring).out
notification = never
request_memory = 2048MB
request_disk = 1024MB
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

rm MakeSFTs.sub datafind.sub test.dag

done
