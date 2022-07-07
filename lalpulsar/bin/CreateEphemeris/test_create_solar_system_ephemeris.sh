# look for files in this directory
export TEMPO2=$PWD
export LALPULSAR_DATADIR=$PWD

# test ephemeris creation against reference file ./earth98.dat in test_create_solar_system_ephemeris.tar.gz
# - Note that this is *not* the same file as lalpulsar/test/earth98.dat (created 2001, commit f79e92da916fb931eeb85ff01bf1510cbb0992f9).
#   lalpulsar_create_solar_system_ephemeris --verbose --test no longer compares equal against this file, presumably because there have been
#   changes to the ephemeris computations since 2001. Instead I have created a new ./earth98.dat with the same parameters as a reference result.
#   Full details are available in the header information for ./earth98.dat
# - The TEMPO2 DE405.1950.2050 ephemeris is taken from https://bitbucket.org/psrsoft/tempo2.git, commit d9133db93d7bac8410cb8445596bf542e43b9317.
lalpulsar_create_solar_system_ephemeris --verbose --test
