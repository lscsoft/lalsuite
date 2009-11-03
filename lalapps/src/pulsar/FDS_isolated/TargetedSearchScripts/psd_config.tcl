

# Configuration Script for the coherent followup pipeline.

foreach {var value} {
	work_dir  "/archive/home/josephb/psd_estimates"

	sft_location_files  {
		H1 "/archive/home/josephb/followScript/TotalListOnNodes_H1"
		L1  "/archive/home/josephb/followScript/TotalListOnNodes_L1"
		}
	
	log_file_dir "/usr1/josephb"

	script_dir  "/archive/home/josephb/PowerFluxFollowUpResamp"
	
	veto_segment_file_name unknown
	segment_file_name unknown
	
	
	} {
	global $var
	set $var [subst $value]
	#puts "$var [set $var]"
	}
