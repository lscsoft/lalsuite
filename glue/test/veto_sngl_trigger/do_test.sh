#!/bin/sh

ligolw_veto_sngl_trigger --trigger-dir triggers --type inspiral --exclude-veto --output-dir results --gps-start-time 800000000 --gps-end-time 800000016 --veto-file veto_definer.xml --ifos H1 --categories 1,3 --dq-dir dq_flags --separate --cumulative

echo 800000003  > wanted
echo 800000004 >> wanted
echo 800000005 >> wanted
echo 800000006 >> wanted
echo 800000007 >> wanted
echo 800000011 >> wanted
echo 800000015 >> wanted

ligolw_print results/H1-CAT1/H1-CAT1-8000/H1-CAT1-800000000-16.xml -t sngl_inspiral -c end_time > result

diff result wanted  >& /dev/null

if [ $? == 1 ]
then
        echo "Cat 1 failed"
        exit 1 
fi

ligolw_print results/H1-CAT1_CUMULATIVE/H1-CAT1_CUMULATIVE-8000/H1-CAT1_CUMULATIVE-800000000-16.xml -t sngl_inspiral -c end_time > result

diff result wanted  >& /dev/null

if [ $? == 1 ]
then
        echo "Cat 1 cumulative failed"
        exit 1 
fi

echo 800000011  > wanted
echo 800000012 >> wanted
echo 800000013 >> wanted

ligolw_print results/H1-CAT3/H1-CAT3-8000/H1-CAT3-800000000-16.xml -t sngl_inspiral -c end_time > result

diff result wanted  >& /dev/null

if [ $? == 1 ]
then
        echo "Cat 3 failed"
        exit 1 
fi

echo 800000011  > wanted

ligolw_print results/H1-CAT3_CUMULATIVE/H1-CAT3_CUMULATIVE-8000/H1-CAT3_CUMULATIVE-800000000-16.xml -t sngl_inspiral -c end_time > result

diff result wanted  >& /dev/null

if [ $? == 1 ]
then
        echo "Cat 3 cumulative failed"
        exit 1 
fi

echo "All tests succeeded."

rm -rf results
rm wanted result

