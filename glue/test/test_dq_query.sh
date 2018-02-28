#!/bin/bash

ligolw_dq_query --segment file://dq_query_testdoc-800000000-16.xml --defined 800000000 | ligolw_print -t segment_definer -c name > result
echo TEST_SEG_0  > wanted
echo TEST_SEG_1 >> wanted
echo TEST_SEG_2 >> wanted
echo TEST_SEG_3 >> wanted

diff result wanted  >& /dev/null

if [ $? == 1 ]
then
        echo "Test 1 failed"
        exit 1 
fi


ligolw_dq_query --segment file://dq_query_testdoc-800000000-16.xml --defined 800000010 | ligolw_print -t segment_definer -c name > result
echo TEST_SEG_2  > wanted
echo TEST_SEG_3 >> wanted

diff result wanted  >& /dev/null

if [ $? == 1 ]
then
        echo "Test 2 failed"
        exit 1 
fi

ligolw_dq_query --segment file://dq_query_testdoc-800000000-16.xml --defined --include-segments H1:TEST_SEG_2 800000010 | ligolw_print -t segment_definer -c name > result
echo TEST_SEG_2  > wanted

diff result wanted  >& /dev/null

if [ $? == 1 ]
then
        echo "Test 3 failed"
        exit 1 
fi

ligolw_dq_query --segment file://dq_query_testdoc-800000000-16.xml --active 800000003 | ligolw_print -t segment_definer -c name > result
echo TEST_SEG_0  > wanted
echo TEST_SEG_3 >> wanted

diff result wanted  >& /dev/null

if [ $? == 1 ]
then
        echo "Test 4 failed"
        exit 1 
fi

ligolw_dq_query --segment file://dq_query_testdoc-800000000-16.xml --active 800000003 --include-segments H1:TEST_SEG_0 | ligolw_print -t segment_definer -c name > result
echo TEST_SEG_0  > wanted

diff result wanted  >& /dev/null

if [ $? == 1 ]
then
        echo "Test 5 failed"
        exit 1 
fi

ligolw_dq_query --segment file://dq_query_testdoc-800000000-16.xml --active 800000003 800000005 | ligolw_print -t segment_definer -c name > result
echo TEST_SEG_0  > wanted
echo TEST_SEG_1 >> wanted
echo TEST_SEG_3 >> wanted

diff result wanted  >& /dev/null

if [ $? == 1 ]
then
        echo "Test 6 failed"
        exit 1 
fi

ligolw_dq_query --segment file://dq_query_testdoc-800000000-16.xml --active --end-pad 5 800000003 | ligolw_print -t segment_definer -c name > result
echo TEST_SEG_0  > wanted
echo TEST_SEG_1 >> wanted
echo TEST_SEG_3 >> wanted

diff result wanted  >& /dev/null

if [ $? == 1 ]
then
        echo "Test 7 failed"
        exit 1 
fi


echo "All tests succeeded"
rm result wanted

