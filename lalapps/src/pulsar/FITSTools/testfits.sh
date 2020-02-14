## create test commands
echo "lalapps_fits_overview 'FITSFileIOTest.fits'"         >  test_lalapps_fits_overview.sh
echo "lalapps_fits_overview 'FITSFileIOTest.fits[array1]'" >> test_lalapps_fits_overview.sh
echo "lalapps_fits_overview 'FITSFileIOTest.fits[array2]'" >> test_lalapps_fits_overview.sh
echo "lalapps_fits_overview 'FITSFileIOTest.fits[array3]'" >> test_lalapps_fits_overview.sh
echo "lalapps_fits_overview 'FITSFileIOTest.fits[array4]'" >> test_lalapps_fits_overview.sh
echo "lalapps_fits_overview 'FITSFileIOTest.fits[table1]'" >> test_lalapps_fits_overview.sh

echo "lalapps_fits_header_list 'FITSFileIOTest.fits'"         >  test_lalapps_fits_header_list.sh
echo "lalapps_fits_header_list 'FITSFileIOTest.fits[array1]'" >> test_lalapps_fits_header_list.sh
echo "lalapps_fits_header_list 'FITSFileIOTest.fits[array2]'" >> test_lalapps_fits_header_list.sh
echo "lalapps_fits_header_list 'FITSFileIOTest.fits[array3]'" >> test_lalapps_fits_header_list.sh
echo "lalapps_fits_header_list 'FITSFileIOTest.fits[array4]'" >> test_lalapps_fits_header_list.sh
echo "lalapps_fits_header_list 'FITSFileIOTest.fits[table1]'" >> test_lalapps_fits_header_list.sh

echo "lalapps_fits_table_list 'FITSFileIOTest.fits[table1]'" > test_lalapps_fits_table_list.sh

echo "lalapps_fits_header_getval 'FITSFileIOTest.fits'         'DATE'"         >  test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits'         'PROGNAME'"     >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits'         'LAL VERSION'"  >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits'         'TESTBOOL'"     >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits'         'TESTUSHRT'"    >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits'         'TESTDBL'"      >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits'         'TESTDBLCMP'"   >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits'         'TESTSTR'"      >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits'         'DATE-OBS'"     >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits'         'DATE-OBS GPS'" >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits[array1]' 'BITPIX'"       >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits[array1]' 'NAXIS1'"       >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits[array2]' 'BITPIX'"       >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits[array2]' 'NAXIS1'"       >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits[array3]' 'BITPIX'"       >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits[array3]' 'NAXIS1'"       >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits[array4]' 'BITPIX'"       >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits[array4]' 'NAXIS1'"       >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits[table1]' 'BITPIX'"       >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits[table1]' 'NAXIS1'"       >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits[table1]' 'TFIELDS'"      >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits[table1]' 'TTYPE1'"       >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits[table1]' 'TFORM1'"       >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits[table1]' 'TTYPE10'"      >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits[table1]' 'TFORM10'"      >> test_lalapps_fits_header_getval.sh
echo "lalapps_fits_header_getval 'FITSFileIOTest.fits[table1]' 'TUNIT10'"      >> test_lalapps_fits_header_getval.sh

## run test commands
for testcmd in *.sh; do
    echo "$SHELL ./$testcmd"
    $SHELL ./$testcmd > stdout-$testcmd.txt 2>stderr.txt
    if [ -s stderr.txt ]; then
        echo "ERROR: $SHELL ./$testcmd should not write to standard error"
        exit 1
    fi
done

## create reference tarball to make updating it a bit easier
mkdir -p newtarball/
cp FITSFileIOTest.fits newtarball/
for file in stdout-*.txt; do
    cp $file newtarball/ref-$file
done
cd newtarball/
tar zcf ../new_testfits.tar.gz *
cd ..
rm -rf newtarball/

## compare standard output to reference results
for file in stdout-*.txt; do
    if ! diff -s $file ref-$file; then
        echo "ERROR: $file and ref-$file should be equal"
        exit 1
    fi
done
