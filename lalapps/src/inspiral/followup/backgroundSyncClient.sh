#!/bin/bash
#
# Sync the ~cbc background directories to here
# Index associates servers to sync paths etc
# 0=LLO, 1=LHO 2=Virgo
localPath=$HOME"/followupbackgrounds/"
tconvert=/ligotools/bin/tconvert
gridsshpath=/opt/ldg-4.7-1/globus/bin/ssh.d/
#Do not edit below this line!
##################################################################
echo $(date)
source $HOME/.bashrc
source $HOME/.ssh/env
export PATH=$gridsshpath:$PATH
declare -a arrayIFOs
arrayIFOs[0]="L1"
arrayIFOs[1]="H1"
arrayIFOs[2]="V1"

declare -a arrayServers
arrayServers[0]="ldas-grid.ligo-la.caltech.edu"
arrayServers[1]="ldas-grid.ligo-wa.caltech.edu"
#arrayServers[2]="ldas-grid.ligo-wa.caltech.edu"

declare -a remoteLog
remoteLog[0]="/archive/home/ctorres/followupbackgrounds/omega/S6/background/latest-run.log"
remoteLog[1]="/archive/home/ctorres/followupbackgrounds/omega/S6/background/latest-run.log"
#remoteLog[2]="/archive/home/ctorres/followupbackgrounds/omega/S6/background/latest-run.log"

#Do NOT put ending slashes on paths we want to sync the
#directory named here to avoid clobbering files
declare -a arraySourcePaths
arraySourcePaths[0]="/archive/home/ctorres/followupbackgrounds/omega"
arraySourcePaths[1]="/archive/home/ctorres/followupbackgrounds/omega"
#arraySourcePaths[2]="/archive/home/cbc/omega/backgrounds"

#Simple error check the length of arraySourcePaths
#should match the length of arrayServers
serverLen=${#arrayServers[*]}

#Loop through all servers check the server log if 
#the jobs have been completed for the day sync them
#otherwise tell the user the background scans are not
#ready to be synced.

#MAIN LOOP
index=0
while [ "$index" -lt "$serverLen" ]
do
    echo "Attempting to sync ${arrayIFOs[$index]} from ${arraySourcePaths[$index]} on ${arrayServers[$index]}"
    #Mkdir if 
    mkdir -p $localPath
    #Sync the information
    echo "Syncing Omega backgrounds last compiled at $endTime from ${arrayServers[$index]}"
    cmd=$(echo rsync -avz  ${arrayServers[$index]}:${arraySourcePaths[$index]} $localPath)
    eval $cmd
    #Check all found log files
    for mylog in $(find $localPath/* | grep latest-run.log); 
    do
	omegaReady=$(tail -2 $mylog | grep "ended at:" &>/dev/null;echo $?)
	endTime=$(tail -1 $mylog | awk -F" " '{print $1}')
	if [ "$omegaReady" -ne "0" ]; then
	    echo "Sync for server may be incomplete! :"${arrayIFOs[$index]}
	    echo "Current GPS Time :$($tconvert now)"
	    echo "The log contents from this server are:"
	    echo "Path to bad log :",$mylog
	    echo $(cat $mylog)
	fi
    done    
    let "index=$index+1"
    echo
done
echo "Follow up background information is available at"
echo $localPath





################
#Scrap Commands#
################
#cmd=$(echo rsync -avz --include=\"${arrayIFOs[$index]}_**\" --exclude=\"*\" ${arrayServers[$index]}:${arraySourcePaths[$index]} $localPath)
