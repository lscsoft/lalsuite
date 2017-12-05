#!/bin/bash
#
# This script is a tempory soltution for the change in the lalapps
#inspiral.c behavior using the hostname to determine where the output
#of this code is written to.  This script will attempt to prevent the
#job from failing and place the output data products back into the
#locations we original designed the followup pipeline for.  Calling it
#without the new directory present creates is.  Calling it again with
#directory there moves the contents of that directory up one level and
#remove the hostname based directory from the disk.  This script will
#also remove the symlinks places in the directory with the dag file.
#
# Config these programs
binHOSTNAME="/bin/hostname"
binMKDIR="/bin/mkdir"
binSED="/bin/sed"
binMV="/bin/mv"
binRM="/bin/rm"
binFIND="/usr/bin/find"
# End config section
#
if [ $# -eq 1 ]; then
    myHostname=`$binHOSTNAME`
    newPath=$1"/"$myHostname"/"
    if [ -d $newPath ]; then
	echo "Moving data products around... Fixing!"
	#Move files up one level
	$binMV $newPath* $1
	$binRM -rf $newPath
	#rm new Symlinks 3 levels higher
	symLinkLevel=$1"/../../"
	find $symLinkLevel -maxdepth 1 -type l -print
	for symLink in $($binFIND $symLinkLevel -maxdepth 1 -type l);do
	    if [ -h $symLink ]; then
		rm -f $symLink
	    fi
	done
    else
	echo "Creating node dependent pathing... Fixing!"
	$binMKDIR --parents $newPath
    fi
else
    echo "Usage: followup_InspiralDataMover.sh ORIGINAL_--output-path_ARG"
fi
exit 0
