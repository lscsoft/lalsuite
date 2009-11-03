#!/bin/bash

# this script will be used as part of 'make test' to make sure that the
# lalapps_heterodyne_pulsar code outputs the correct result. It will run 
# the code in all its 4 modes and compare the results to standard
# archived files 

CODENAME=./lalapps_heterodyne_pulsar

FRAMEFILE=framefiles.tar.gz
DATASTART=875206560
DATAEND=`expr $DATASTART + 24 \* 60 \* 60`
DETECTOR=H1
CHANNEL=LSC-DARM_ERR
FKNEE=0.25

# sample rates (frames are 1024Hz)
SRATE1=1024
SRATE2=1
SRATE3=1/60

# create a pulsar par file
PSRNAME=J0000+0000
FREQ=245.678910
FDOT=-9.87654321e-12
RA=00:00:00.0
DEC=00:00:00.0
PEPOCH=53966.22281462963
PFILE=$PSRNAME.par

if [ -f $PFILE ]; then
	rm -f $PFILE
fi

echo PSR    $PSRNAME > $PFILE
echo F0     $FREQ >> $PFILE
echo F1     $FDOT >> $PFILE
echo RAJ    $RA >> $PFILE
echo DECJ   $DEC >> $PFILE
echo PEPOCH $PEPOCH >> $PFILE

if [ $? != "0" ]; then
	echo Error writing parameter file!
	exit 2
fi

# create slightly offset parameter file (for testing update mode)
FREQOFF=245.679
FDOTOFF=-9.9e-12
PEPOCHOFF=53966.0
PFILEOFF=${PSRNAME}_offset.par

if [ -f $PFILEOFF ]; then
        rm -f $PFILEOFF
fi

echo PSR    $PSRNAME > $PFILEOFF
echo F0     $FREQOFF >> $PFILEOFF
echo F1     $FDOTOFF >> $PFILEOFF
echo RAJ    $RA >> $PFILEOFF
echo DECJ   $DEC >> $PFILEOFF
echo PEPOCH $PEPOCH >> $PFILEOFF

if [ $? != "0" ]; then
        echo Error writing parameter file!
        exit 2
fi

# set ephemeris file
EEPHEM=$LAL_PREFIX/share/lal/earth05-09.dat
if [ ! -f $EEPHEM ]; then
	echo Error! Earth ephemeris file does not exist!
	exit 2
fi

SEPHEM=$LAL_PREFIX/share/lal/sun05-09.dat
if [ ! -f $SEPHEM ]; then
	echo Error! Sun ephemeris file does not exist!
	exit 2
fi 

# get current location
LOCATION=`pwd`
if [ $? != "0" ]; then
	echo Error! Could not set the current path!
	exit 2
fi

# create a frame cache file in the format used by the code

# check file doesn't already exit
if [ -f cachefile ]; then
	rm -f cachefile
fi

# make directory to contain the frames and unpack the tar file
# first check if it already exists
if [ -d ${LOCATION}/framedir ]; then
	rm -f ${LOCATION}/framedir/*
	rmdir ${LOCATION}/framedir
fi

mkdir ${LOCATION}/framedir
if [ $? != "0" ]; then
	echo Error. Could not create frame directory
	exit 2
fi

if [ ! -f $FRAMEFILE ]; then
	echo Error. Frame tar file does not exist!
	exit 2
fi

tar -xvzf $FRAMEFILE
if [ $? != "0" ]; then
	echo Error. Could not unpack the frame files!
	exit 2
fi
FILELIST=`tar -tzf $FRAMEFILE`
mv $FILELIST ${LOCATION}/framedir 

# use make_frame_cache to make a frame cache file
if [ ! -f make_frame_cache ]; then
	echo Error! make_frame_cache does not exist!
	exit 2
fi

./make_frame_cache --frame-dir ${LOCATION}/framedir --gps-start-time $DATASTART --gps-end-time $DATAEND --output-file cachefile
if [ $? != "0" ]; then
	echo Could not create the cache file!
	exit 2
fi

# create segment file

# check file doesn't already exit
if [ -f segfile ]; then
	rm -f segfile
fi

# make 24 segments of two minutes length each spanning a day
for (( i=0; i<24; i++ )); do
	SEGNUM=`expr $i + 1`
	SEGSTART=`expr $DATASTART + $i \* 3600`
	SEGEND=`expr $SEGSTART + 120`
	echo $SEGNUM $SEGSTART $SEGEND 120 >> segfile
done
if [ $? != "0" ]; then
	echo Could not create the segment file!
	exit 2
fi

# make output directory in format that the code like
OUTDIR=$LOCATION/$DATASTART-$DATAEND

# check if it exists first
if [ -d $OUTDIR ]; then
	rm -f ${OUTDIR}/*
	rmdir $OUTDIR
fi

mkdir $OUTDIR
if [ $? != "0" ]; then
	echo Could not create the output directory
	exit 2
fi

#################### COARSE HETERODYNES ########################

# run code in coarse heterodyne mode (outputing to a text file)
echo Performing coarse heterodyne - mode 0 - and outputting to text file 
$CODENAME --heterodyne-flag 0 --ifo $DETECTOR --pulsar $PSRNAME --param-file $PFILE --sample-rate $SRATE1 --resample-rate $SRATE2 --filter-knee $FKNEE --data-file $LOCATION/cachefile --seg-file $LOCATION/segfile --channel $CHANNEL --output-dir $OUTDIR --freq-factor 2

# check the exit status of the code
if [ $? != "0" ]; then
	echo lalapps_heterodyne_pulsar exited with error $?!
	exit 2
fi

# check that the expected file got output
COARSEFILE=$OUTDIR/coarsehet_${PSRNAME}_${DETECTOR}_${DATASTART}-${DATAEND}
if [ ! -f $COARSEFILE ]; then
	echo Error! Code has not output a coarse heterodyne file
	exit 2
fi

# move file so that we can create a binary file of the same name
mv $COARSEFILE $COARSEFILE.txt

# run code in coarse heterodyne mode again (outputing to a binary file)
echo Performing coarse heterodyne - mode 0 - and outputting to binary file
$CODENAME --heterodyne-flag 0 --ifo $DETECTOR --pulsar $PSRNAME --param-file $PFILE --sample-rate $SRATE1 --resample-rate $SRATE2 --filter-knee $FKNEE --data-file $LOCATION/cachefile --seg-file $LOCATION/segfile --channel $CHANNEL --output-dir $OUTDIR --binary-output --freq-factor 2

# check the exit status of the code
if [ $? != "0" ]; then
        echo lalapps_heterodyne_pulsar exited with error $?!
        exit 2
fi

# check that the expected file got output
if [ ! -f $COARSEFILE ]; then
        echo Error! Code has not output a coarse heterodyne file
        exit 2
fi

mv $COARSEFILE $COARSEFILE.bin 

# run code in coarse heterodyne mode again, but this time with the offset par file
echo Performing coarse heterodyne - mode 0 - with offset parameter file
$CODENAME --heterodyne-flag 0 --ifo $DETECTOR --pulsar $PSRNAME --param-file $PFILEOFF --sample-rate $SRATE1 --resample-rate $SRATE2 --filter-knee $FKNEE --data-file $LOCATION/cachefile --seg-file $LOCATION/segfile --channel $CHANNEL --output-dir $OUTDIR --freq-factor 2

# check the exit status of the code
if [ $? != "0" ]; then
        echo lalapps_heterodyne_pulsar exited with error $?!
        exit 2
fi

# check that the expected file got output
if [ ! -f $COARSEFILE ]; then
        echo Error! Code has not output a coarse heterodyne file
        exit 2
fi

mv $COARSEFILE $COARSEFILE.off

# set calibration files
RESPFILE=$LOCATION/H1response.txt

################### FINE HETERODYNES #######################

# now perform the fine heterodyne (first using the txt file)
echo Performing fine heterodyne - mode 1 - using text file
$CODENAME --ephem-earth-file $EEPHEM --ephem-sun-file $SEPHEM --heterodyne-flag 1 --ifo $DETECTOR --pulsar $PSRNAME --param-file $PFILE --sample-rate $SRATE2 --resample-rate $SRATE3 --filter-knee $FKNEE --data-file $COARSEFILE.txt --output-dir $OUTDIR --channel $CHANNEL --seg-file $LOCATION/segfile --freq-factor 2 --calibrate --response-file $RESPFILE --stddev-thresh 5

# check the exit status of the code
if [ $? != "0" ]; then
        echo lalapps_heterodyne_pulsar exited with error $?!
        exit 2
fi

# check that it produced the right file
FINEFILE=$OUTDIR/finehet_${PSRNAME}_${DETECTOR}
if [ ! -f $FINEFILE ]; then
	echo Error! Code has not output a fine heterodyned file
	exit 2
fi

# move file
mv $FINEFILE $FINEFILE.txt

# now perform the fine heterodyne (first using the binary file)
echo Performing fine heterodyne - mode 1 - using binary file
$CODENAME --ephem-earth-file $EEPHEM --ephem-sun-file $SEPHEM --heterodyne-flag 1 --ifo $DETECTOR --pulsar $PSRNAME --param-file $PFILE --sample-rate $SRATE2 --resample-rate $SRATE3 --filter-knee $FKNEE --data-file $COARSEFILE.bin --binary-input --output-dir $OUTDIR --channel $CHANNEL --seg-file $LOCATION/segfile --freq-factor 2 --calibrate --response-file $RESPFILE --stddev-thresh 5

# check the exit status of the code
if [ $? != "0" ]; then
        echo lalapps_heterodyne_pulsar exited with error $?!
        exit 2
fi

# check that it produced the right file
if [ ! -f $FINEFILE ]; then
        echo Error! Code has not output a fine heterodyned file
        exit 2
fi

# move file
mv $FINEFILE $FINEFILE.bin

# now perform the fine heterodyne with the updating that with offset parameter file
echo Performing fine heterodyne - mode 2 - using update from offset parameter file
$CODENAME --ephem-earth-file $EEPHEM --ephem-sun-file $SEPHEM --heterodyne-flag 2 --ifo $DETECTOR --pulsar $PSRNAME --param-file $PFILEOFF --param-file-update $PFILE --sample-rate $SRATE2 --resample-rate $SRATE3 --filter-knee $FKNEE --data-file $COARSEFILE.off --output-dir $OUTDIR --channel $CHANNEL --seg-file $LOCATION/segfile --freq-factor 2 --calibrate --response-file $RESPFILE --stddev-thresh 5

# check the exit status of the code
if [ $? != "0" ]; then
        echo lalapps_heterodyne_pulsar exited with error $?!
        exit 2
fi

# check that it produced the right file
if [ ! -f $FINEFILE ]; then
        echo Error! Code has not output a fine heterodyned file
        exit 2
fi

# move file
mv $FINEFILE $FINEFILE.off

# now perform the fine heterodyne with the offset parameter file (no update)
echo Performing fine heterodyne - mode 1 - using offset parameter file
$CODENAME --ephem-earth-file $EEPHEM --ephem-sun-file $SEPHEM --heterodyne-flag 1 --ifo $DETECTOR --pulsar $PSRNAME --param-file $PFILEOFF --sample-rate $SRATE2 --resample-rate $SRATE3 --filter-knee $FKNEE --data-file $COARSEFILE.off --output-dir $OUTDIR --channel $CHANNEL --seg-file $LOCATION/segfile --freq-factor 2 --calibrate --response-file $RESPFILE --stddev-thresh 5

# check the exit status of the code
if [ $? != "0" ]; then
        echo lalapps_heterodyne_pulsar exited with error $?!
        exit 2
fi

# check that it produced the right file
if [ ! -f $FINEFILE ]; then
        echo Error! Code has not output a fine heterodyned file
        exit 2
fi

# move file
mv $FINEFILE $FINEFILE.off2

################### HETERODYNE ALL IN ONE #############
# now perform the heterodyne in one go (mode 3)
echo Performing entire heterodyne in one go - mode 3
$CODENAME --ephem-earth-file $EEPHEM --ephem-sun-file $SEPHEM --heterodyne-flag 3 --ifo $DETECTOR --pulsar $PSRNAME --param-file $PFILE --sample-rate $SRATE1 --resample-rate $SRATE3 --filter-knee $FKNEE --data-file $LOCATION/cachefile --output-dir $OUTDIR --channel $CHANNEL --seg-file $LOCATION/segfile --freq-factor 2 --calibrate --response-file $RESPFILE --stddev-thresh 5

# check the exit status of the code
if [ $? != "0" ]; then
	echo lalapps_heterodyne_pulsar exited with error $?!
	exit 2
fi

# check that it produced the right file
if [ ! -f $FINEFILE ]; then
	echo Error! Cde has not output a fine heterodyned file
	exit 2
fi

# move file
mv $FINEFILE $FINEFILE.full

################### REHETERODYNE THE ALREADY FINE HETERODYNED FILE #####
echo Performing updating heterodyne of already fine heterodyned data
$CODENAME --ephem-earth-file $EEPHEM --ephem-sun-file $SEPHEM --heterodyne-flag 4 --ifo $DETECTOR --pulsar $PSRNAME --param-file $PFILEOFF --param-file-update $PFILE --sample-rate $SRATE3 --resample-rate $SRATE3 --filter-knee 0 --data-file $FINEFILE.off2 --output-dir $OUTDIR --channel $CHANNEL --seg-file $LOCATION/segfile --freq-factor 2 --stddev-thresh 5

# check the exit status of the code
if [ $? != "0" ]; then
        echo lalapps_heterodyne_pulsar exited with error $?!
        exit 2
fi

# check that it produced the right file
if [ ! -f $FINEFILE ]; then
        echo Error! Code has not output a fine heterodyned file
        exit 2
fi

###### CHECK THAT ALL VERSION PRODUCE MATCH THE PRE-MADE REFERENCES #####  
echo Comparing outputs with reference files

# check reference files exist
if [ ! -f simulated_reference.txt ]; then
	echo Error. Reference file does not exist!
	exit 2
fi

if [ ! -f finehet_J0000+0000_H1 ]; then
	echo Error. Reference file does not exist!
        exit 2
fi

./compare_test_outputs

if [ $? != 0 ]; then
	echo Error. Files do not match references!
	exit 2
fi

################### CLEAN UP ##########################
echo Cleaning up directory.

# remove cache file
rm -f cachefile

# remove segment file
rm -f segfile

# remove parameter files
rm -f $PFILE
rm -f $PFILEOFF

# remove upacked frame files
rm -f ${LOCATION}/framedir/*
rmdir ${LOCATION}/framedir

# remove files produced during heterodyne
rm -f ${OUTDIR}/*
rmdir $OUTDIR

if [ $? != "0" ]; then
	echo Error. Something went wrong during clean up!
	exit 2
fi

# exit with all being well :)
echo All is well with the world :D

exit 0
