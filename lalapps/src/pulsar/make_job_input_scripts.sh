#!/bin/bash

#a1=S2-H2.60
#a2=L 
#a3=L1
#a4=S2-H2.60-CAL-lists
#a5=10

a1=$1
a2=$2
a3=$3
a4=$4
a5=$5

mkdir $a4
njobs=$a5
jobtotal=`awk '{ N+=$1 ; print N }' $a1 | tail -1`
#echo will make $jobtotal segments

echo "#"\!"/bin/bash"
jobn=0
while [ $jobn -lt $njobs ]
do
  cat $a1 \
   | awk '{ N+=$1 ; print $1, $2, $3, N}' \
   | awk '{ if ($4 > '$jobtotal'*'$jobn'/'$njobs' ) if ($4 <= '$jobtotal'*('$jobn'+1)/'$njobs') print $1, $2, $3, $4, '$jobn'}' \
   | awk '{printf "%s" "%05d" "%s\n", "LSCdataFind --server=dataserver.phys.uwm.edu --observatory '$a2' --type '$a3'_RDS_C02_LX --gps-start-time "$2" --gps-end-time "$2+$1*$3" --url-type file --match localhost | sed s,file://localhost,, >> '$a4'/jobdata.", ""$5"", ".ffl" ; for (sft=0; sft < $1; ++sft) printf "%s" "%05d""%s\n",  "echo "$2+sft*$3" "$3" >> '$a4'/jobtimes.", ""$5"", "" }'
  let "jobn+=1" 
done
echo "exit 0"

exit 0
             
