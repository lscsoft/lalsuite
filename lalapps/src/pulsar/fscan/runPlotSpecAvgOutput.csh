#!/bin/tcsh
#
# run plotSpecAvgOutput(inputFileName,outputFileName,chanName,tStart,tEnd,fStart,fEnd,effTBase,deltaFTicks,medBins)
# The path to plotSpecAvgOutput must be in the MATLABPATH  env variable.
# See plotSpecAvgOutput.m for meaning of command line arguments.
$MATLAB_ROOT/bin/matlab -nodisplay -nodesktop -nojvm << EOF
plotSpecAvgOutput('$1','$2','$3','$4','$5','$6','$7','$8','$9','$10')
exit
EOF

