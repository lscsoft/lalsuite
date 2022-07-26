# look for files in this directory
export TEMPO2=$PWD

# test ephemeris creation against reference files ./ref_tdb_2000-2019.dat.gz and ./ref_te405_2000-2019.dat.gz
# - These are copied from lalpulsar/lib (created 2012, commit 50e428e00ad970617d6f302cb87a36fd8fb62cb4) and uncompressed.
# - The TDB.1950.2050 and TIMEEPH_short.te405 ephemerides are taken from https://bitbucket.org/psrsoft/tempo2.git, commit d9133db93d7bac8410cb8445596bf542e43b9317.
lalapps_create_time_correction_ephemeris --verbose --ephem-type TDB --output-path . --start 630720000 --end 1262736000 --interval 14400
lalapps_create_time_correction_ephemeris --verbose --ephem-type TCB --output-path . --start 630720000 --end 1262736000 --interval 14400

# compare ephemeris files
for file in tdb_2000-2019.dat te405_2000-2019.dat; do
    echo
    echo "Comparing ${file} and ref_${file}"
    cat "${file}" | sed '/^#/d' > new.txt
    cat "ref_${file}" | sed '/^#/d' > ref.txt
    paste new.txt ref.txt | awk -v "rmsmax=1e-12" '
        BEGINFILE { rms = 0 }
        FNR == 1 { if ( $1 != $5 || $2 != $6 || $3 != $7 || $4 != $8 ) { print "Header mismatch:", $0; exit ( 1 ) } }
        FNR > 1 { rms += ($1 - $2)^2 }
        ENDFILE { rms = sqrt(rms / (NR - 1)); printf "rms=%g rmsmax=%g\n", rms, rmsmax; exit ( rms <= rmsmax ? 0 : 1 ) }
    '
done
