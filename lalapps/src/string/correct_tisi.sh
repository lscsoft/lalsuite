#!/bin/sh

if [ $# != 1 ]; then 
    echo "correct_tisi.sh ERROR : the time slide file should be given in argument"
    exit 1
fi

echo "This macro will remove the time slides where the difference"
echo "of time offset between V1 and L1 appears multiple times"
echo "Please wait it can take some time..."

tisi_file_gz="$1.xml.gz"

# check the tisi file
if [ ! -e ${tisi_file_gz} ]; then
    echo "correct_tisi.sh ERROR : ${tisi_file_gz} does not exist"
    exit 1
fi

gunzip ${tisi_file_gz}
tisi_file="$1.xml"

new_tisi_file="new_${tisi_file}"
tmp_tisi_file="tmp_${tisi_file}"
rm -f ${new_tisi_file} ${tmp_tisi_file} 

# read the file line by line
config=0; number_of_timeslide=0; previous_offset=9999999;
while read line; do
    
    # re-write time slides based on L1 offsets
    if echo ${line} | grep -q 'time_slide:time_slide_id:'; then
	if echo ${line} | grep -q '"L1"'; then
	    L1offset=`echo ${line} | cut -d"," -f4`
	    V1offset=`echo "0-($L1offset)" | bc`

	    if [ ! ${L1offset} == ${previous_offset} ]; then
		sed '/<Stream Name="time_slide:table" Type="Local" Delimiter=",">/a\ "H1","time_slide:time_slide_id:'$config'","process:process_id:0",0,' ${new_tisi_file} > ${tmp_tisi_file}
		mv ${tmp_tisi_file} ${new_tisi_file}
		sed '/<Stream Name="time_slide:table" Type="Local" Delimiter=",">/a\ "H2","time_slide:time_slide_id:'$config'","process:process_id:0",0,' ${new_tisi_file} > ${tmp_tisi_file}
		mv ${tmp_tisi_file} ${new_tisi_file}
		sed '/<Stream Name="time_slide:table" Type="Local" Delimiter=",">/a\ "L1","time_slide:time_slide_id:'$config'","process:process_id:0",'$L1offset',' ${new_tisi_file} > ${tmp_tisi_file}
		mv ${tmp_tisi_file} ${new_tisi_file}
		sed '/<Stream Name="time_slide:table" Type="Local" Delimiter=",">/a\ "V1","time_slide:time_slide_id:'$config'","process:process_id:0",'$V1offset',' ${new_tisi_file} > ${tmp_tisi_file}
		mv ${tmp_tisi_file} ${new_tisi_file}
		echo "ti-si ${config}: L1 offset=$L1offset, V1 offset=$V1offset"
		
		let 'config+=1'
		previous_offset=${L1offset}
	    fi
	fi
	let 'number_of_timeslide+=1'
    else
	echo ${line} >> ${new_tisi_file}
    fi
done < ${tisi_file}

let 'number_of_timeslide/=4'

echo "The original ${number_of_timeslide} time-slides have been reduced to ${config} independent time-slides"


mv ${new_tisi_file} ${tisi_file}
gzip ${tisi_file}

exit 0
