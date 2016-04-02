#!/bin/tcsh
# The path to plotSpecAvgOutput must be in the MATLABPATH env variable.
# See plotSpecAvgOutput.m for meaning of command line arguments.
if ($#argv != 10) then
    $MATLAB_ROOT/bin/matlab -nodisplay -nodesktop -nojvm << EOF
help plotSpecAvgOutput
exit
EOF
else
    $MATLAB_ROOT/bin/matlab -nodisplay -nodesktop -nojvm << EOF
plotSpecAvgOutput('$1','$2','$3','$4','$5','$6','$7','$8','$9','$10')
exit
EOF
endif
echo
