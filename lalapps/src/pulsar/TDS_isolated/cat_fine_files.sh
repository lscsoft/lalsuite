#!/bin/sh

# Matt Pitkin 10/05/07


# cat_fine_files.sh - this is a bash script that will take in any number of directories as command
# line arguments, list the files in the directories, and concatenate any that have the same name
# and output them to the last directory give in the command line arguments (which it will create if
# it doesn't already exist)

# the main point of this script is to concatenate the fine heterodyned data files output in the
# known pulsar search - although it can be a bit more versatile

# get the number of command line inputs
numdirs=$#

# check that their are more than three arguments (two input directories and one output directory)
if [ $numdirs -lt 3 ] # -lt means less than
then
  echo "Their must be three or more input arguments"
  exit 1
fi

count=1

# check that all bar the last directory actually exist (if the last one doesn't then make it)
for i in $@
do
  if [ ! -d $i ] && [ $count -lt $numdirs ] # -d means does the directory exist
  then
    echo "Directory $i does not exist, boo!"
  else
    if [ $count -eq $numdirs ]
    then
      # get name of output directory for later
      outputdir=$i
    fi
    
    if [ ! -d $i ] && [ $count -eq $numdirs ] # -a means AND
    then
      # create final directory
      echo "Creating output directory $i"
      mkdir $i
    fi
  fi

  count=`expr $count + 1` # remember that you need spaces in the expression
done

count=1

# go though input directories, test whether their are files with the same name, and conatenate to
# the output directory if this is true
for i in $@
do
  # list files in first directory
  if [ $count -eq 1 ]
  then
    filelist=`ls $i`
    count=`expr $count + 1`
    
    # output the files in the first dir to the output file
    for a in $filelist
    do
      outputfile=$outputdir/$a
      
      cat $i/$a >> $outputfile # >> means append
    done
    
    continue
  fi
  
  # output files in the rest of the directories
  if [ $count -gt 1 ] && [ $count -lt $numdirs ]
  then
    filelisttwo=`ls $i`
    
    for a in $filelist
    do
      for b in $filelisttwo
      do
        # check if files have the same name
        if [ "$a" = "$b" ]
        then
          # concatenate files out to output file
          outputfile=$outputdir/$a
          cat $i/$b >> $outputfile 
        fi
      done
    done
  fi
  
  count=`expr $count + 1`
done

exit 0
