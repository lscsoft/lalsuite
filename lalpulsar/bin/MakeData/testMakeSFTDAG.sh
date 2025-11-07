## test exit status
status=0

## common variables
Tsft=1800
seg1_tstart=1257741529
seg1_tend=`echo "${seg1_tstart} + ${Tsft}" | bc`
seg1_sft1="${seg1_tstart}"
seg2_tstart=1257743330
seg2_tend=`echo "${seg2_tstart} + 1.5*${Tsft}" | bc | xargs printf '%0.0f'`
seg2_sft1_half_overlap="${seg2_tstart}"
seg2_sft2_half_overlap=`echo "${seg2_tstart} + 0.5*${Tsft}" | bc | xargs printf '%0.0f'`
seg2_sft1_no_overlap=`echo "${seg2_tstart} + 0.25*${Tsft}" | bc | xargs printf '%0.0f'`
frame_dur=4096
frame2_tstart=`echo "${seg1_tstart} + 4096" | bc`
fmin=10
Band=1990
chan1="H1:GDS-CALIB_STRAIN_CLEAN"
chan2="H1:GDS-CALIB_STRAIN"
acctgtag="ligo.sim.o4.cw.explore.test"
acctgusr="albert.einstein"
MSFTpath="/tmp/path/to"
cachepath="${PWD}/frames.cache"

## channels names in SFT names
chan1sft=`echo "${chan1}" | sed 's/^H1://;s/[-_]//g'`
chan2sft=`echo "${chan2}" | sed 's/^H1://;s/[-_]//g'`

## common segment file
segs="${PWD}/segs"
echo "${seg1_tstart} ${seg1_tend}" > ${segs}
echo "${seg2_tstart} ${seg2_tend}" >> ${segs}

## common cache file
echo -e "H\tH1_HOFT_C00\t${seg1_tstart}\t${frame_dur}\tH-H1_HOFT_C00-${seg1_tstart}-${frame_dur}.gwf" > ${cachepath}
echo -e "H\tH1_HOFT_C00\t${frame2_tstart}\t${frame_dur}\tH-H1_HOFT_C00-${frame2_tstart}-${frame_dur}.gwf" >> ${cachepath}

## variable for to silence gwdatafind error
export GWDATAFIND_SERVER=TEST

## common unit test function
unittest() {

    # first argument is name of unit test
    # remaining arguments are passed to lalpulsar_MakeSFTDAG
    local name="$1"
    shift

    echo
    echo "############### unit test: ${name} ###############"
    echo

    # run lalpulsar_MakeSFTDAG in unit test directory
    local testdir="./${name}_unittest"
    rm -rf ${testdir}
    mkdir ${testdir}
    cd ${testdir}
    echo "Running lalpulsar_MakeSFTDAG..."
    ../lalpulsar_MakeSFTDAG "$@"
    output_code=$?
    echo "Output code of lalpulsar_MakeSFTDAG $@ is: $output_code"
    if [[ $output_code -ne 0 ]]; then
        if [[ $output_code -ne 77 ]]; then
            echo "ERROR: something failed when running 'lalpulsar_MakeSFTDAG $@'"
            exit 1
        else
            echo "Skipping test - no gwpy installed"
            return 77
        fi
    fi
    cd ..
    echo

    # store contents of all files in unit test directory
    local output="./${name}_output.txt"
    find ${testdir} -type f | LC_ALL=C sort | xargs grep -H . > ${output}

    # compare output to reference
    local output_ref="./${name}_output_ref.txt"
    touch ${output_ref}
    if ! diff ${output_ref} ${output} >/dev/null 2>&1; then
        echo "ERROR: output of unit test ${name} differs from reference as follows:"
        echo "--------------------------------------------------------------------------------"
        diff ${output_ref} ${output} || true
        echo "--------------------------------------------------------------------------------"
        echo
        status=1

        # create new reference tarball
        cp -f ${output} ${output_ref}
        rm -f testMakeSFTDAG.tar.gz
        tar zcf testMakeSFTDAG.tar.gz *_output_ref.txt

    else
        echo "INFO: output of unit test ${name} is identical to reference:"
        echo "--------------------------------------------------------------------------------"
        cat ${output}
        echo "--------------------------------------------------------------------------------"
        echo
    fi

}

## function to grep test output for strings
greptest() {

    # first argument is name of unit test
    # remaining arguents are strings to grep
    local name="$1"
    shift

    # grep test output for strings
    local output="./${name}_output.txt"
    while test "X$1" != X; do
        if grep -q -e "$1" "${output}"; then
            echo "INFO: output of unit test ${name} contains the string '$1'"
        else
            echo "ERROR: output of unit test ${name} should contain the string '$1'"
            status=1
        fi
        echo
        shift
    done

}

############### unit tests ###############

## public SFTs
unittest public_SFTs \
    -O 4 -K DEV -R 1 \
    -f test.dag -G TEST -d H1_HOFT_C00 -k 7 -T ${Tsft} -p . \
    -N ${chan1} -F ${fmin} -B ${Band} -w hann -P 0.5 -m 1 \
    -A ${acctgtag} -U ${acctgusr} \
    -g ${segs} -J ${MSFTpath} --movesfts-path=${MSFTpath} \
    -e ${cachepath}
if [[ $? -ne 77 ]]; then
    greptest public_SFTs \
         "-O 4 -K DEV -R 1" \
         "-w hann -P 0.5"
fi

## two SFTs per job
unittest two_SFTs_per_job \
    -O 4 -K DEV -R 1 \
    -f test.dag -G TEST -d H1_HOFT_C00 -k 7 -T ${Tsft} -p . \
    -N ${chan1} -F ${fmin} -B ${Band} -w hann -P 0.5 -m 2 \
    -A ${acctgtag} -U ${acctgusr} \
    -g ${segs} -J ${MSFTpath} --movesfts-path=${MSFTpath} \
    -e ${cachepath}
if [[ $? -ne 77 ]]; then
    greptest two_SFTs_per_job \
         "-O 4 -K DEV -R 1" \
         "-w hann -P 0.5"
fi

## two channels
unittest two_channels \
    -O 4 -K DEV -R 1 \
    -f test.dag -G TEST -d H1_HOFT_C00 -k 7 -T ${Tsft} -p dir1 dir2 \
    -N ${chan1} ${chan2} -F ${fmin} -B ${Band} -w hann -P 0.5 -m 1 \
    -A ${acctgtag} -U ${acctgusr} \
    -g ${segs} -J ${MSFTpath} --movesfts-path=${MSFTpath} \
    -e ${cachepath}
if [[ $? -ne 77 ]]; then
    greptest two_channels \
         "-O 4 -K DEV -R 1" \
         "-w hann -P 0.5"
fi

## private SFTs
unittest private_SFTs \
    -O 0 -X private \
    -f test.dag -G TEST -d H1_HOFT_C00 -k 7 -T ${Tsft} -p . \
    -N ${chan1} -F ${fmin} -B ${Band} -w hann -P 0.5 -m 1 \
    -A ${acctgtag} -U ${acctgusr} \
    -g ${segs} -J ${MSFTpath} --movesfts-path=${MSFTpath} \
    -e ${cachepath}
if [[ $? -ne 77 ]]; then
    greptest private_SFTs \
        "-O 0 -X private" \
        "-w hann -P 0.5"
fi

## default window
unittest default_window \
    -O 4 -K DEV -R 1 \
    -f test.dag -G TEST -d H1_HOFT_C00 -k 7 -T ${Tsft} -p . \
    -N ${chan1} -F ${fmin} -B ${Band} -m 1 \
    -A ${acctgtag} -U ${acctgusr} \
    -g ${segs} -J ${MSFTpath} --movesfts-path=${MSFTpath} \
    -e ${cachepath}
if [[ $? -ne 77 ]]; then
    greptest default_window \
        "-O 4 -K DEV -R 1" \
        "-w tukey -r 0.001"
fi

## Tukey window with parameter 0.001
unittest Tukey_window \
    -O 4 -K DEV -R 1 \
    -f test.dag -G TEST -d H1_HOFT_C00 -k 7 -T ${Tsft} -p . \
    -N ${chan1} -F ${fmin} -B ${Band} -w tukey:0.001 -m 1 \
    -A ${acctgtag} -U ${acctgusr} \
    -g ${segs} -J ${MSFTpath} --movesfts-path=${MSFTpath} \
    -e ${cachepath}
if [[ $? -ne 77 ]]; then
    greptest Tukey_window \
        "-O 4 -K DEV -R 1" \
        "-w tukey -r 0.001"
fi

## Tukey window with parameter 0.5
unittest Tukey_window_2 \
    -O 4 -K DEV -R 1 \
    -f test.dag -G TEST -d H1_HOFT_C00 -k 7 -T ${Tsft} -p . \
    -N ${chan1} -F ${fmin} -B ${Band} -w tukey:0.5 -m 1 \
    -A ${acctgtag} -U ${acctgusr} \
    -g ${segs} -J ${MSFTpath} --movesfts-path=${MSFTpath} \
    -e ${cachepath}
if [[ $? -ne 77 ]]; then
    greptest Tukey_window_2 \
         "-O 4 -K DEV -R 1" \
         "-w tukey -r 0.5"
fi

if test -f testMakeSFTDAG.tar.gz; then
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "INFO: new reference result tarball has been generated: '${PWD}/testMakeSFTDAG.tar.gz'"
    echo "INFO: copy to the source directory of 'testMakeSFTDAG.sh' to update the test"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo
fi

exit ${status}
