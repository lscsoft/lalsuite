#!/bin/sh
# simple (naive?) script to produce the H+ and H- frames

if [ $# != 2 ]; then 
    exit 1
fi

gps_start=$1
gps_end=$2
outdir=/archive/home/frobinet/jw1string/H1H2-coherent

mkdir -p ${outdir}

# H1 and H2 cache files
ligo_data_find --observatory H --url-type file --gps-start-time ${gps_start} --gps-end-time ${gps_end}  --lal-cache --type H1_RDS_C03_L2 > ${outdir}/H1.cache
ligo_data_find --observatory H --url-type file --gps-start-time ${gps_start} --gps-end-time ${gps_end}  --lal-cache --type H2_RDS_C03_L2 > ${outdir}/H2.cache

lalapps_StringAddFrame ${gps_start} ${gps_end} ${outdir}

exit 0
