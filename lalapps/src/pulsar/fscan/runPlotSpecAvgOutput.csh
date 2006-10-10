#!/bin/tcsh
#
#
# run plotSpecAvgOutput(inputFileName,outputFileName,chanName,tStart,tEnd,fStart,fEnd,effTBase,deltaFTicks)
$MATLAB_ROOT/bin/matlab -nodisplay -nodesktop -nojvm << EOF
plotSpecAvgOutput('$1','$2','$3','$4','$5','$6','$7','$8','$9')
exit
EOF

