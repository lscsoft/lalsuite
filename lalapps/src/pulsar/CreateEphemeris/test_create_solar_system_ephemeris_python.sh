# this test requires the 'jplephem' module
if ${PYTHON} -c "import jplephem"; then
    echo "jplephem is available; running test"
else
    echo "jplephem is not available; skipping test"
    exit 77
fi

# test ephemeris creation against reference file ./earth98.dat
# - For details of ./earth98.dat see ./test_create_solar_system_ephemeris.sh
# - The binary SPK file ./de405-1998.bsp was creating using the command:
#      python3 -m jplephem excerpt 1997/12/15 1999/1/15 https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp de405-1998.bsp
lalpulsar_create_solar_system_ephemeris_python --ephemeris ./de405-1998.bsp --output-file ./earth98-python.dat --gps-start 566784000 --num-years 1 --interval 4 --target EARTH

# strip comments from ephemeris files, put both in same (newer) format of one GPS time/position/velocity/acceleration entry per line
cat earth98.dat | sed '/^#/d' | sed '1d' | awk '{ ORS=""; print $1, $2, $3; getline; print " ", $1, $2, $3; getline; print " ", $1, $2, $3; getline; print " ", $1, "\n"}' | awk 'NR <= 2192 { print }' > earth98.txt
cat earth98-python.dat | sed '/^#/d' | sed '1d' > earth98-python.txt

# compare ephemeris files
echo
echo "Comparing earth98.txt and earth98-python.txt"
for i in `seq 1 10`; do
    case $i in
        1)      rmsmax='0';;            # GPS time
        2|3|4)  rmsmax='1e-8';;         # position
        5|6|7)  rmsmax='1e-15';;        # velocity
        8|9|10) rmsmax='1e-20';;        # acceleration
    esac
    paste earth98.txt earth98-python.txt | awk -v "i=$i" -v "rmsmax=$rmsmax" '
        BEGINFILE { j = i + 10; rms = 0 }
        { rms += ($i - $j)^2 }
        ENDFILE { rms = sqrt(rms / NR); printf "i=%i rms=%g rmsmax=%g\n", i, rms, rmsmax; exit ( rms <= rmsmax ? 0 : 1 ) }
    '
done
