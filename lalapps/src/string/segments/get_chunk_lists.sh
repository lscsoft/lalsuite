#!/bin/bash


c=1

while read line; do
    chunk_start=`echo $line | awk '{print $2}'`
    chunk_end=`echo $line | awk '{print $3}'`
    if [ "$chunk_start" = "" ]; then continue; fi
    echo "Building segment lists for chunk ${c}..."
    echo $chunk_start $chunk_end > ./tmp
    segexpr 'intersection(tmp,h1.seg)' -include index,duration > ./S5-H1segments-cat1_${c}.txt
    segexpr 'intersection(tmp,h2.seg)' -include index,duration > ./S5-H2segments-cat1_${c}.txt
    segexpr 'intersection(tmp,l1.seg)' -include index,duration > ./S5-L1segments-cat1_${c}.txt
    segexpr 'intersection(tmp,v1.seg)' -include index,duration > ./S5-V1segments-cat1_${c}.txt
        
    h1_lt=`segsum ./S5-H1segments-cat1_${c}.txt`
    h2_lt=`segsum ./S5-H2segments-cat1_${c}.txt`
    l1_lt=`segsum ./S5-L1segments-cat1_${c}.txt`
    v1_lt=`segsum ./S5-V1segments-cat1_${c}.txt`
    echo "H1 livetime=${h1_lt} sec"
    echo "H2 livetime=${h2_lt} sec"
    echo "L1 livetime=${l1_lt} sec"
    echo "V1 livetime=${v1_lt} sec"

    let "c+=1"
done < chunks.txt

rm -f ./tmp

exit 0
