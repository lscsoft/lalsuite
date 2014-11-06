#!/bin/bash

# Test the gracedb command line client
#

: ${GRACEDB:="gracedb"}
: ${TEST_DATA_DIR:=$(dirname $0)/data}
: ${GRACEDB_SERVICE_URL:="https://moe.phys.uwm.edu/branson/api/"}

export GRACEDB_SERVICE_URL

N=0
NSUCC=0
NFAIL=0
NERR=0

function recordTest() {
    TESTNAME=$1
    RETCODE=$2
    OUT="$3"
    case $RETCODE in
        0)  NSUCC=$[$NSUCC+1]
            MESSAGE="Succeeded"
            ;;
        1)  NFAIL=$[$NFAIL+1]
            MESSAGE="Failed $OUT"
            ;;
        *)  NERR=$[$NERR+1]
            MESSAGE="Error $OUT"
            ;;
    esac
    N=$[$N+1]
    echo $TESTNAME $MESSAGE
}

function showStats() {
    echo "Success: " $NSUCC
    echo "Fail:    " $NFAIL
    echo "Error:   " $NERR
    echo "Total:   " $N
}

OUT="$(${GRACEDB} ping 2>&1)"
recordTest "ping" "$?" "$OUT"

# Try creating events of various types
#

OUT="$(${GRACEDB} Test gstlal LowMass $TEST_DATA_DIR/cbc-lm.xml 2>&1)"
RETCODE=$?
recordTest "create gstlal" "$RETCODE" "$OUT"

# Remember the LM event.  We will use it later
#
if [ $RETCODE -eq 0 ]
then
    GRACEID=$OUT
else
    GRACEID="NOID"
fi

OUTFILE=$(mktemp /tmp/tmp.XXXXXXXXX)
${GRACEDB} Test MBTAOnline $TEST_DATA_DIR/cbc-mbta.gwf >$OUTFILE 2>&1
recordTest "create MBTA" "$?" "$(cat $OUTFILE)"
rm $OUTFILE

OUTFILE=$(mktemp /tmp/tmp.XXXXXXXXX)
${GRACEDB} Test CWB $TEST_DATA_DIR/burst-cwb.txt >$OUTFILE 2>&1
recordTest "create CWB"  "$?" "$(cat $OUTFILE)"
rm $OUTFILE


# Try a simple search
#
OUTFILE=$(mktemp /tmp/tmp.XXXXXXXXX)
${GRACEDB} search $GRACEID >$OUTFILE 2>&1
recordTest "search $GRACEID" "$?" "$(cat $OUTFILE)"
rm $OUTFILE

# Try a simple search with ligolw output.
#
OUTFILE=$(mktemp /tmp/tmp.XXXXXXXXX)
${GRACEDB} search $GRACEID --ligolw >$OUTFILE 2>&1
recordTest "search $GRACEID --ligolw" "$?" "$(cat $OUTFILE)"
rm $OUTFILE

# Make sure FAR of created LM event is correct.
#
OUTFILE=$(mktemp /tmp/tmp.XXXXXXXXX)
${GRACEDB} search "--columns=gpstime" $GRACEID > $OUTFILE 2>&1
RETCODE=$?
if [ $RETCODE == 0 ]
then
    if [ "$(grep -v '#' <$OUTFILE)" == 971609248 ]
    then
        RETCODE=0
    else
        RETCODE=1
    fi
fi
recordTest "verify GPS time $GRACEID" "$RETCODE" "$(cat $OUTFILE)"
rm $OUTFILE


# Replace LM event with new data
#
OUTFILE=$(mktemp /tmp/tmp.XXXXXXXXX)
${GRACEDB} replace $GRACEID $TEST_DATA_DIR/cbc-lm2.xml >$OUTFILE 2>&1
RETCODE=$?
recordTest "replace $GRACEID" "$RETCODE" "$(cat $OUTFILE)"
rm $OUTFILE


# Make sure FAR of replaced LM event is correct.
#
OUTFILE=$(mktemp /tmp/tmp.XXXXXXXXX)
${GRACEDB} search "--columns=gpstime" $GRACEID > $OUTFILE 2>&1
RETCODE=$?
if [ $RETCODE == 0 ]
then
    if [ "$(grep -v '#' <$OUTFILE)" == 971609249 ]
    then
        RETCODE=0
    else
        RETCODE=1
    fi
fi
recordTest "verify new GPS time $GRACEID" "$RETCODE" "$(cat $OUTFILE)"
rm $OUTFILE

# Upload a file
#
OUTFILE=$(mktemp /tmp/tmp.XXXXXXXXX)
${GRACEDB} --tag-name=tag_test upload $GRACEID "$TEST_DATA_DIR/upload.data.gz" > $OUTFILE 2>&1
recordTest "upload file $GRACEID" "$?" "$(cat $OUTFILE)"
rm $OUTFILE

# Download that uploaded file
#
DOWNLOAD=$(mktemp /tmp/tmp.XXXXXXXXX)
${GRACEDB} download $GRACEID "upload.data.gz" - > $DOWNLOAD 2>&1
recordTest "download file" "$?" "$(cat $DOWNLOAD)"

# Verify that the uploaded file and downloaded file were the same
#
cmp --silent "$DOWNLOAD" "$TEST_DATA_DIR/upload.data.gz"
recordTest "verify uploaded file" "$?" "$(cat $DOWNLOAD)"
rm $DOWNLOAD

# Log
#
OUTFILE=$(mktemp /tmp/tmp.XXXXXXXXX)
${GRACEDB} log $GRACEID "test the" "logging" > $OUTFILE 2>&1
recordTest "log message" "$?" "$(cat $OUTFILE)"
rm $OUTFILE

# Label
#
OUTFILE=$(mktemp /tmp/tmp.XXXXXXXXX)
${GRACEDB} label $GRACEID DQV > $OUTFILE 2>&1
recordTest "label $GRACEID" "$?" "$(cat $OUTFILE)"
rm $OUTFILE

# Verify labelling
#
OUTFILE=$(mktemp /tmp/tmp.XXXXXXXXX)
${GRACEDB} search $GRACEID |grep DQV > $OUTFILE 2>&1
recordTest "verify label" "$?" "$(cat $OUTFILE)"
rm $OUTFILE

# Tag
#
OUTFILE=$(mktemp /tmp/tmp.XXXXXXXXX)
${GRACEDB} tag $GRACEID 1 tag_test > $OUTFILE 2>&1
recordTest "tag event" "$?" "$(cat $OUTFILE)"
rm $OUTFILE

# Delete tag
#
OUTFILE=$(mktemp /tmp/tmp.XXXXXXXXX)
${GRACEDB} delete_tag $GRACEID 1 tag_test > $OUTFILE 2>&1
recordTest "delete tag" "$?" "$(cat $OUTFILE)"
rm $OUTFILE

showStats

