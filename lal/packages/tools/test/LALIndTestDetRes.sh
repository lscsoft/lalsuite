#!/bin/sh
# Shell script that runs independent test of F_+ and F_x using LALIndependentTestDetResponse
echo "./LALIndependentTestDetResponse -c indTestDetResLHO.cfg"
./LALIndependentTestDetResponse -c indTestDetResLHO.cfg || exit
echo "./LALIndependentTestDetResponse -c indTestDetResLLO.cfg"
./LALIndependentTestDetResponse -c indTestDetResLLO.cfg || exit
echo "./LALIndependentTestDetResponse -c indTestDetResGEO.cfg"
./LALIndependentTestDetResponse -c indTestDetResGEO.cfg || exit
exit
