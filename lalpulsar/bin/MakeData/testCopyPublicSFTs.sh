echo
echo "===== Test directory construction on dummy SFTs ====="
echo

## make dummy source SFTs
mkdir -p ./src
for src_SFT in `cat src_SFT_files.txt`; do
    src_subdir=`echo ${src_SFT} | cut -d- -f 3 | cut -c 9-`
    mkdir -p "./src/${src_subdir}"
    echo "dummy" > "./src/${src_subdir}/${src_SFT}"
done

## copy SFTs
mkdir -p ./dest
lalpulsar_CopyPublicSFTs --no-validate --processes 2 ./src ./dest

## check SFTs were copied to the right place
find dest -name '*.sft' | sort > ./dest_SFT_files.txt
diff -s ./ref_dest_SFT_files.txt ./dest_SFT_files.txt

echo
echo "===== Test validation of SFTs ====="
echo

## make good/bad SFTs
SFTwrite

## rename good SFTs - not filename is assumed consitent with contents
mkdir -p ./src-good
for i in 1 2 3; do
    src_SFT="H-1_H1_1800SFT_O4DEV+V1+Choft+WRECT-100020030${i}-1800.sft"
    cp "./SFT-test${i}" "./src-good/${src_SFT}"
done

## copy good SFTs
mkdir -p ./dest-good
lalpulsar_CopyPublicSFTs ./src-good ./dest-good

## check good SFTs were copied to the right place
for i in 1 2 3; do
    src_SFT="H-1_H1_1800SFT_O4DEV+V1+Choft+WRECT-100020030${i}-1800.sft"
    cmp "./src-good/${src_SFT}" "./dest-good/H1_1800SFT_O4DEV+V1+Choft+WRECT_BROADBAND/GPS1000M/${src_SFT}"
    echo "OK: ./src-good/${src_SFT} ->  ./dest-good/H1_1800SFT_O4DEV+V1+Choft+WRECT_BROADBAND/GPS1000M/${src_SFT}"
done

## rename a bad SFT
mkdir -p ./src-bad
cp ./SFT-bad1 ./src-bad/H-1_H1_1800SFT_O4DEV+V1+Cbadt+WRECT-999999999-1800.sft

## copy bad SFT - should fail
mkdir -p ./dest-bad
if lalpulsar_CopyPublicSFTs ./src-bad ./dest-bad; then
    echo "ERROR: lalpulsar_CopyPublicSFTs ./src-bad ./dest-bad should have failed"
fi
