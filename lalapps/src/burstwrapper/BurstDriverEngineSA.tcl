############################################################################
## The function injFileCopy is used to make a copy of injection data file
## with the right substitutions for different IFOs
############################################################################

proc injFileCopy { data start_time IFO2 } {

    if { [ catch { exec cp ${data}.ilwd.$start_time ${data}${IFO2}.ilwd } err ] } {
	puts "Error: can't copy file!"
	puts $err
	exit -1
    }

    set fin [ open ${data}.ilwd.$start_time "r" ]
    set fout [ open ${data}${IFO2}.ilwd "w" ]

    gets $fin line

    while { [ eof $fin ] == 0 } {
	regsub -all $data $line $data$IFO2 oline
	puts $fout $oline

	gets $fin line
    }

    close $fin
    close $fout
}

############################################################################
## Modifies channel name to get ADC data
############################################################################
proc mkchannel { chan stime } {
    regsub -all -- {:} $chan {\\:} chalias
    set chalias ${chalias}::AdcData:$stime:0:Frame
}

############################################################################
## Replaces : with \: in chan
############################################################################
proc mkchannelnt { chan } {
    regsub -all -- {:} $chan {\:} chalias
    return $chalias
}

############################################################################
## Get data from URL
## if outFile is present, saves data to file; otherwise, returns data
############################################################################
proc dataDownload { URL { outFile -1 } } {

	if { [ regexp {http://([^/]+)/.*} $URL junk IP ] } {
		
	if { [ regexp {http://[^/]+(/.*)} $URL junk loc ] } {
	
		set sid [ socket $IP 80 ]
		fconfigure $sid -translation binary
		puts $sid "GET $loc\r\n"
	        flush $sid

		if { $outFile != -1 } {
		    set fid [open $outFile w]
		    fconfigure $fid -translation binary
		    fcopy $sid $fid
		    close $fid
		    close $sid
		} else {
		    set out [read $sid]
		    close $sid
		    return $out
		}

	} else {
	error "Can't find result location in URL $URL"	
	}
	} else {
	 error "Can't find IP address in URL $URL"
	}

}

############################################################################
############################################################################
############################################################################
############################################################################

############################################################################
## Initialize variables
############################################################################
set injAmpTmp $injAmp
set LOCALFILE 0
set FastCacheDone 0

set EndSegTime 0
set cjid 0
set CondorFiles ""

## Make sure sampling frequency is an integer
set samplingF [ expr int($samplingF) ]

## If channel is a comma separated list, split into list named channels
set channels [ split $channel "," ]

## Create list of IFOs from channels list, using first letter of every channel
set IFOs [ list ]
foreach channel $channels {
    lappend IFOs [ string range $channel 0 1 ]
}



############################################################################
## Handles case where PrebinFile starts with LOCAL:
############################################################################
if { [ regexp {LOCAL:(.*)} $PrebinFile junk ZeFile ] } {

    ## Hardcoded number of nodes in the cluster
    set NNODES 210

    ## Remove LOCAL: prefix
    set PrebinFile $ZeFile

    ## Create file on local computer
    if { [ file exists $PrebinFile ] == 0 } {
	file mkdir $PrebinFile
    }

    ## Check that path starts with /usr1
    ## Also save path without .usr1/ prefix in PrebinFileSUFFIX
    if { [ regexp {/usr1/(.*)} $ZeFile junk PrebinFileSUFFIX ] == 0 } {
	error "PrebinFile doesn't start with /usr1/"
    }

    ## Loop over nodes to create path
    for { set ni 1 } { $ni < $NNODES } { incr ni 1 } {

	## Check that node exists
	if { [ file exists "/data/node${ni}" ] } {

	    ## Check that path doesn't exist
	    if { [ file exists "/data/node${ni}/$PrebinFileSUFFIX" ] == 0 } {

		## Create path
		puts "Creating /data/node${ni}/$PrebinFileSUFFIX"
		file mkdir "/data/node${ni}/$PrebinFileSUFFIX"

		## Make path writable universally
		if { [ catch { exec chmod 777 "/data/node${ni}/$PrebinFileSUFFIX" } err ] } {
		    puts "WARNING: can't change mode of /data/node${ni}/$PrebinFileSUFFIX to 777; LDAS might not be able to write there!!"
		}

	    } else {
		## path already exist; dangerous for data corruption
		puts "WARNING: /data/node${ni}/$PrebinFileSUFFIX already exists!"
	    }
	}

    ## end for loop:
    }

} else {

    ########################################################################
    ## PrebinFile isn't LOCAL
    ########################################################################

    ## Create path on local computer
    if { [ file exists $PrebinFile ] == 0 } {
	puts "Creating $PrebinFile"
	file mkdir $PrebinFile
    } else {
	puts "WARNING: $PrebinFile already exists!"
    }

    ## Make universally writable
    if { [ catch { exec chmod 777 $PrebinFile } err ] } {
	puts "WARNING: can't change mode of $PrebinFile to 777; LDAS might not be able to write there!!"
    }

}

########################################################################
## Determine the kind of output
########################################################################
set fileOutput 0
set doSplit 0

switch -exact -- $OutputType {
    0 { set binputput "" }
    1 { set binoutput "binOutput," ; set doSplit 1 }
    2 { set fileOutput 1 }
    default { error "Invalue OutputType: $OutputType" }
}

########################################################################
## burstoutput is added before first argument of filterParams
## if $BurstOutput == 1, noLALBurstOutput is passed, and the
## parameter estimation function is not used.
########################################################################
if { $BurstOutput == 0 } {
    set burstoutput ""
} else {
    set burstoutput "noLALBurstOutput,"
}


########################################################################
## If the code is called with 2 arguments, 2nd argument is used
## as the blockId where to start processing
########################################################################
if { [ llength $argv ] == 2 } {
    set doblock 1
    set blockId [ lindex $argv 1 ]
} else {
    set doblock 0
    set blockId 0
}

########################################################################
## Input data
########################################################################

## eval to handle variable substitution
eval set dataType $dataType

## switch dataType
switch -exact -- $dataType {

    0 { 
	## case 0: get data from playgndURL
	set playgndURL "http://www.ligo.caltech.edu/~jsylvest/S2/coin/Playgnd.txt" 
	set SegmentList [ dataDownload $playgndURL ]
	set SegmentList [ split $SegmentList "\n" ]
    }

    1 { error "Production data mode not yet supported" }

    default { 
	## default is used to pass filenames which contains
	## the start times of all segments

	set NSL 0

	## In case we have a comma separated list of filenames, split
	set SegmentFileList [ split $dataType "," ]

	## foreach file:
	foreach file $SegmentFileList {

	    ## read the times
	    set SegmentListTmp [SegRead $file]

	    if { $NSL < 1 } {
		## if first file to be read
		set SegmentList $SegmentListTmp
	    } else {
		## other files are intersected with first file
		set SegmentList [SegIntersection $SegmentList $SegmentListTmp]
	    }
	    
	    incr NSL 1
	}

        ## Cut to play ground
        if { [ info exists PlaygroundOnly ] } {
            set SegmentList [ SegPlaygroundMask $SegmentList ]
        }

    }
}


##########################################################################
##########################################################################
##########################################################################
## Loop over segments
##########################################################################
##########################################################################
##########################################################################

set jobRetry 0

## index of segment being process
set bId 0 

set do1 1

## begin loop
foreach l1 $SegmentList {

    if { $dataType <= 1 } {

	## If got data from URL, reformat
	if { [ regexp {([0-9\.]+)\s+([0-9\.]+)\s+([01])\s+([01])\s+([01]).*} $l1 junk start_time end_time H1 H2 L1 ] == 0 } {
	    continue
	}

    } else {

	## set start and end times
	set start_time [lindex $l1 0]
	set end_time   [lindex $l1 1]  

	## add one second to end time
	incr end_time  1

	## set I to 1, for I in IFOs list
	foreach I $IFOs {
		set $I 1
	}
    }

    ## make sure start_time is integer
    set start_time [ expr int($start_time) ]

    ## increment segment index
    incr bId 1

    ## run loop to check we have data from all IFOs
    set gotIFO 0
    foreach I $IFOs {
	set gotIFO [ set $I ]
	if { $gotIFO == 0 } {
	    ## data is missing
	    break
	}
    }


    ## Check that we have data (gotIFO), and that block index is
    ## after first requested block ($bId >= $blockId)
    if { $gotIFO && ( $doblock == 0 || $bId >= $blockId ) } {

	## begin loop over all segments of size duration within segment
	while { [ expr $start_time + $duration ] <= $end_time } {

	    ## get random number
	    set seed_return [ expr srand($start_time) ]

	    ## test if we want to skip some segments
	    if { $do1 == $do1inN } {
		set do1 0
	    }
	    incr do1 1

	    ## we want to process this segment
	    if { $do1 == 1 } {

		set channelNumber 0

		## loop over all channels
		foreach channel $channels {

		    incr channelNumber 1

		    set IFO [ string index $channel 0 ]
		    set IFO2 [ string range $channel 0 1 ]
		    set ifo2 [ string tolower $IFO2 ]

		    ## get segment Id from Peter's SegID function
		    set PbId [SegID S2 $IFO2 $start_time]

		    puts "Channel $channel: processing $start_time for $duration seconds (block $bId - Peter's block Id $PbId)"


		    
		    ##########################################################
		    ## Parse prefilters
		    ##########################################################
		    if { [ string length $prefilters ] == 0 } {
			## no prefilters
			set dcfilters ""
			set responseFile ""
		    } else {

			## split comma separated list of prefilters
			set flist [ split $prefilters "," ]
			if { [ llength $flist ] == 0 } {
			    error "Invalid filter list: $prefilters"
			}

			## string with response files to get
			set responseFile ""

			## string with filters actions for datacondAPI
			set dcfilters ""

			## filter index in datacond
			set fcount 1

			## flag
			set got1 0

			## check for "a" coefficients of first filter in list
			set fa [ lindex $flist 0 ]
			append fa "_a.ilwd"
			if { [ file exists $fa ] } {
			    ## local file available, add to responseFile
			    append responseFile "$fa push fa$fcount\n"
			    incr got1 1
			} else {
			    ## no local file, add LDAS filtroot prefix
			    append responseFile "${filtroot}/$fa push fa$fcount\n"
			    incr got1 1
			}

			## check for "b" coefficients, first filter in list
			set fb [ lindex $flist 0 ]
			append fb "_b.ilwd"
			if { [ file exists $fb ] } {
			    append responseFile "$fb push fb$fcount\n"
			    incr got1 1
			} else {
			    append responseFile "${filtroot}/$fb push fb$fcount\n"
			    incr got1 1
			}

			## raise an error is a or b coefficients are missing
			if { $got1 != 2 } {
			    error "Can't file $fa or $fb"
			}

			## add filtering action to datacondAPI algorithms
			append dcfilters "\t\tgwchn = linfilt(fb$fcount,fa$fcount,gwchn);\n"

			## move to other filters
			incr fcount 1

			if { [ llength $flist ] > 1 } {

			    ## loop over filters
			    foreach f [ lrange $flist 1 end ] {

				set got1 0

				## do a coefficients
				set fa $f
				append fa "_a.ilwd"
				if { [ file exists $fa ] } {
				    append responseFile "$fa push fa$fcount\n"
				    incr got1 1
				} else {
				    append responseFile "${filtroot}/$fa push fa$fcount\n"
				    incr got1 1
				}

				## do b coefficients
				set fb $f
				append fb "_b.ilwd"
				if { [ file exists $fb ] } {
				    append responseFile "$fb push fb$fcount\n"
				    incr got1 1
				} else {
				    append responseFile "${filtroot}/$fb push fb$fcount\n"
				    incr got1 1
				}

				## error is a or b is missing
				if { $got1 != 2 } {
				    error "Can't file $fa or $fb"
				}

				## datacond action
				append dcfilters "\t\tgwchn = linfilt(fb$fcount,fa$fcount,gwchn);\n"

				incr fcount 1

			    }

			}
		    }


		    ####################################################################
		    ## Parse prefiltersB
		    ####################################################################
		    if { [ string length $prefiltersB ] > 0 } {

			## sblit comma separated list
			set flist [ split $prefiltersB "," ]

			set got1 0
			
			# look at first filter
			set fb [ lindex $flist 0 ]

			## substitute variable names
			eval set fb $fb 

			## check for filtername.ilwd
			append fb ".ilwd"
			if { [ file exists $fb ] } {
			    ## have local file
			    append responseFile "$fb push fb$fcount\n"
			    incr got1 1
			} else {
			    ## get file from LDAS
			    append responseFile "${filtroot}/$fb push fb$fcount\n"
			    incr got1 1
			}

			## datacond action
			append dcfilters "\t\tgwchn = linfilt(fb$fcount,gwchn);\n"

			## move to next filter
			incr fcount 1

			if { [ llength $flist ] > 1 } {

			    ## loop over filters
			    foreach f [ lrange $flist 1 end ] {
				set got1 0
				set fb $f

				## replace variable names
				eval set fb $fb 
	    
				## loop for file fb.ilwd
				append fb ".ilwd"
				if { [ file exists $fb ] } {
				    ## local file
				    append responseFile "$fb push fb$fcount\n"
				    incr got1 1
				} else {
				    ## need LDAS file
				    append responseFile "${filtroot}/$fb push fb$fcount\n"
				    incr got1 1
				}

				## datacond action
				append dcfilters "\t\tgwchn = linfilt(fb$fcount,gwchn);\n"

				incr fcount 1
			    }
			} 
		    }

		    ## replace gwchn in dcfilters with hpf, save in hpfilters
		    regsub -all "gwchn" $dcfilters "hpf" hpfilters

		    ## replace gwchn in dcfilters with hcf, save in hcfilters
		    regsub -all "gwchn" $dcfilters "hcf" hcfilters


		    ################################################################################
		    ## Parse waveforms
		    ################################################################################
		    if { [ string length $waveforms ] == 0 } {
			## no waveforms
			set injWave ""
			set dcwave ""
		    } else {

			## split comma separated list
			set wlist [ split $waveforms "," ]

			if { [ llength $wlist ] == 0 } {
			    error "Invalid waveform list: $waveforms"
			}

			set wcount 1
			set got0 0
			set got1 0
			set got2 0

			## list of waveforms
			set injWave ""

			## datacond actions for waveforms
			set dcwave ""
    
			## look at 1st waveform
			set ff [ lindex $wlist 0 ]

			## do plus polarization
			set fplus [ lindex $wlist 0 ]

			## look for ilwd file
			append fplus "_p.ilwd"
			if { [ file exists $fplus ] } {
			    ## got local file
			    append responseFile "$fplus push ihp$wcount\n"
			    set got1 1
			}

			## do cross polarization
			set fcross [ lindex $wlist 0 ]

			## look for ilwd file
			append fcross "_c.ilwd"
			if { [ file exists $fcross ] } {
			    ## local file
			    append responseFile "$fcross push ihc$wcount\n"
			    set got2 1
			}

			## have not plus and no cross polarizations
			if { $got1 == 0 && $got2 == 0 } {
			    error "Can't file $fplus or $fcross"
			}

			## one polarization is missing, need file with Zeros
			if { $got0 == 0 && $got1 == 0 } {
			    set f "Zeros.ilwd"
			    if { [ file exists $f ] } {
				append responseFile "Zeros.ilwd push Zero\n"
				append dcwave "\t\tzero = float(Zero)"
				set got0 1
			    }
			}

			if { $got0 == 0 && $got2 == 0 } {
			    set f "Zeros.ilwd"
			    if { [ file exists $f ] } {
				append responseFile "Zeros.ilwd push Zero\n"
				append dcwave "\t\tzero = float(Zero)"
				set got0 1
			    }
			}


			if { $got1 != 0 } {

			    ## We have a plus polarization:
			    ## add datacond code
			    append dcwave "\t\tihpt = tseries(ihp$wcount,$samplingF.0,$start_time);\n"
			    append dcwave "\t\thpf = respfilt(ihpt,h1resp,h1gain,h1cavfacf,h1oloopf);\n"
			    append dcwave $hpfilters
			    append dcwave "\t\thpfr = float(hpf);\n"
			    set st1 "\t\toutput(hpfr,_,_,$ff"
			    append dcwave [ append st1 "_p,$ff" "_p);\n" ]
			} else {

			    ## no plus polarization
			    ## add zeros
			    set st1 "\t\toutput(zero,_,_,$ff"
			    append dcwave [ append st1 "_p,$ff" "_p);\n" ]
			}

			if { $got2 != 0 } {
			    ## We have a cross polarization:
			    ## add datacond code
			    append dcwave "\t\tihct = tseries(ihc$wcount,$samplingF.0,$start_time);\n"
			    append dcwave "\t\thcf = respfilt(ihct,h1resp,h1gain,h1cavfacf,h1oloopf);\n"
			    append dcwave $hcfilters
			    append dcwave "\t\thcfr = float(hcf);\n"
			    set st1 "\t\toutput(hcfr,_,_,$ff"
			    append dcwave [ append st1 "_c,$ff" "_c);\n" ]
			} else {

			    ## no cross polarization
			    ## add zeros
			    set st1 "\t\toutput(zero,_,_,$ff"
			    append dcwave [ append st1 "_c,$ff" "_c);\n" ]
			}

			incr wcount 1
    

			set w _[ lindex $wlist 0 ]
			set injWave [ append w "_" ]

			## do the other waveforms
			if { [ llength $wlist ] > 1 } {

			    ## use ( to start list of waveforms
			    set injWave "($injWave"

			    ## loop over waveforms
			    foreach w [ lrange $wlist 1 end ] {

				append injWave ",_" $w "_"

				set got1 0
				set got2 0

				## plus polarization
				set fplus $w
				append fplus "_p.ilwd"
				if { [ file exists $fplus ] } {
				    append responseFile "$fplus push ihp$wcount\n"
				    set got1 1
				}

				## cross polarization
				set fcross $w
				append fcross "_c.ilwd"
				if { [ file exists $fcross ] } {
				    append responseFile "$fcross push ihc$wcount\n"
				    set got2 1
				}

				## error if no polarization
				if { $got1 == 0 && $got2 == 0 } {
				    error "Can't file $fplus or $fcross"
				}

				## handle missing polarization with zeros
				if { $got0 == 0 && $got1 == 0 } {
				    set f "Zeros.ilwd"
				    if { [ file exists $f ] } {
					append responseFile "Zeros.ilwd push Zero\n"
					append dcwave "\t\tzero = float(Zero)"
					set got0 1
				    }
				}

				if { $got0 == 0 && $got2 == 0 } {
				    set f "Zeros.ilwd"
				    if { [ file exists $f ] } {
					append responseFile "Zeros.ilwd push Zero\n"
					append dcwave "\t\tzero = float(Zero)"
					set got0 1
				    }
				}

				
				if { $got1 != 0 } {
				    ## plus polarization datacond actions
				    append dcwave "\t\tihpt = tseries(ihp$wcount,$samplingF.0,$start_time);\n"
				    append dcwave "\t\thpf = respfilt(ihpt,h1resp,h1gain,h1cavfacf,h1oloopf);\n"
				    append dcwave $hpfilters
				    append dcwave "\t\thpfr = float(hpf);\n"
				    set st1 "\t\toutput(hpfr,_,_,$w"
				    append dcwave [ append st1 "_p,$w" "_p);\n" ]
				} else {
				    ## zeros
				    set st1 "\t\toutput(zero,_,_,$w"
				    append dcwave [ append st1 "_p,$w" "_p);\n" ]
				}

				if { $got2 != 0 } {
				    ## cross polarization datacond actions
				    append dcwave "\t\tihct = tseries(ihc$wcount,$samplingF.0,$start_time);\n"
				    append dcwave "\t\thcf = respfilt(ihct,h1resp,h1gain,h1cavfacf,h1oloopf);\n"
				    append dcwave $hcfilters
				    append dcwave "\t\thcfr = float(hcf);\n"
				    set st1 "\t\toutput(hcfr,_,_,$w"
				    append dcwave [ append st1 "_c,$w" "_c);\n" ]
				} else {
				    ## zeros
				    set st1 "\t\toutput(zero,_,_,$w"
				    append dcwave [ append st1 "_c,$w" "_c);\n" ]
				}
				
				incr wcount 1
				
			    }

			    ## close list with )
			    append injWave ")"

			}
		    }


		    ################################################################################
		    ## Parse injType for injections from spatial distributions
		    ################################################################################
		    if { [ info exists injType ] == 1 } {
			switch -exact -- $injType {
			    1 { 
				
				## Do Shell Injection

				## run ShellInjection1 to produce injection data
				if { $channelNumber == 1 } {
				    if { [ catch { exec ShellInjection1 $injAmpTmp $injNtot $injN $duration } cout ] } {
					puts $cout
					exit -1
				    }
				}

				## set the variables to output of ShellInjection1
				set injAmp "__hamp__"
				set injAlpha "__alpha__"
				set injDelta "__delta__"
				set injPsi "__psi__"
				set injTimes "__injTime__"

				## save a backup copy of injected stuff
				if { $injSave == 1 } {
				    if { [ catch { exec cp hamp.ilwd hamp.ilwd.$start_time } cout ] } {
					puts $cout
					exit -1
				    }
				    if { [ catch { exec cp alpha.ilwd alpha.ilwd.$start_time } cout ] } {
					puts $cout
					exit -1
				    }
				    if { [ catch { exec cp delta.ilwd delta.ilwd.$start_time } cout ] } {
					puts $cout
					exit -1
				    }
				    if { [ catch { exec cp psi.ilwd psi.ilwd.$start_time } cout ] } {
					puts $cout
					exit -1
				    }
				    if { [ catch { exec cp injTime.ilwd injTime.ilwd.$start_time } cout ] } {
					puts $cout
					exit -1
				    }

				}
			    }


			    3 { 

				## run ShellInjection
				if { $channelNumber == 1 } {
				    if { [ catch { exec ShellInjection $injAmpTmp $injNtot $injN $duration } cout ] } {
					puts $cout
					exit -1
				    }
				}

				## set variables to output of ShellInjection
				set injAmp "__hamp__"
				set injAlpha "__alpha__"
				set injDelta "__delta__"
				set injPsi "__psi__"
				set injTimes "__injTime__"

				## Save backup copy in injected stuff
				if { $injSave == 1 } {
				    if { [ catch { exec cp hamp.ilwd hamp.ilwd.$start_time } cout ] } {
					puts $cout
					exit -1
				    }
				    if { [ catch { exec cp alpha.ilwd alpha.ilwd.$start_time } cout ] } {
					puts $cout
					exit -1
				    }
				    if { [ catch { exec cp delta.ilwd delta.ilwd.$start_time } cout ] } {
					puts $cout
					exit -1
				    }
				    if { [ catch { exec cp psi.ilwd psi.ilwd.$start_time } cout ] } {
					puts $cout
					exit -1
				    }
				    if { [ catch { exec cp injTime.ilwd injTime.ilwd.$start_time } cout ] } {
					puts $cout
					exit -1
				    }
				    
				}
			    }


			    2 { 

				## Use files generated by another driver

				set skip 0

				## Loop until file becomes available
				set now0 [ clock seconds ]
				while { [ file exists hamp.ilwd.$start_time ] == 0 } {
				    set now [ clock seconds ]
				    if { [ expr $now - $now0 ] > 1000 } {
					puts "Couldn't get injection files; skipping"
					set skip 1
					break
				    }
				    after 10000
				}

				if { $skip } {
				    continue
				}

				## Point to the right data to use
				set injAmp "__hamp${IFO2}__"
				set injAlpha "__alpha${IFO2}__"
				set injDelta "__delta${IFO2}__"
				set injPsi "__psi${IFO2}__"
				set injTimes "__injTime${IFO2}__"

				## copy files from other driver
				injFileCopy "hamp" $start_time ${IFO2} 	    
				injFileCopy "alpha" $start_time ${IFO2} 	    
				injFileCopy "delta" $start_time ${IFO2} 	    
				injFileCopy "psi" $start_time ${IFO2} 	    
				injFileCopy "injTime" $start_time ${IFO2} 	    

			    }

			    default {
				puts "invalid injType: $injType"
				exit -1
			    }
			}
		    }


		    ################################################################################
		    ## Handle case where injection data is a matrix (__FILE__ notation)
		    ################################################################################
		    if { [ regexp {__([^_]+)__} $injAmp junk file ] } {
			append responseFile ${file}".ilwd pass\n"
		    }
		    if { [ regexp {__([^_]+)__} $injAlpha junk file ] } {
			append responseFile ${file}".ilwd pass\n"
		    }
		    if { [ regexp {__([^_]+)__} $injDelta junk file ] } {
			append responseFile ${file}".ilwd pass\n"
		    }
		    if { [ regexp {__([^_]+)__} $injPsi junk file ] } {
			append responseFile ${file}".ilwd pass\n"
		    }
		    if { [ regexp {__([^_]+)__} $injTimes junk file ] } {
			append responseFile ${file}".ilwd pass\n"
		    }


		    ################################################################################
		    ## Parse ETG parameters
		    ################################################################################

		    if { [ info exists ETGParameters ] == 0 } { 
			## Explicit list of parameters
			if { $ETG == "TFCLUSTERS" } {
			    set TFCchn [ set channel ]
			    set ETGParams "$TFCchn,$TFCThr,$TFCWin,$TFCThrMethod,$TFCdata,$TFCTRez,$TFCFmin,$TFCFmax,$TFCalpha,$TFCsigma,$TFCdelta"
			}
		    } else {
			## ETGParameters passed as a matrix file (__FILE__ notation)

			## eval to do variable substitution
			eval set ETGParameters $ETGParameters

			if { [ regexp {__(.+)__} $ETGParameters junk file ] } {
			    ## add file to responseFile
			    append responseFile "${file}.ilwd pass\n"
			}

			set ETGParams $ETGParameters
		    }




		    ################################################################################
		    ## Prepare job variables
		    ################################################################################

		    ## output file
		    set binFile "$PrebinFile/job$IFO2.$userName.$start_time.$bId.bin"

		    set binoutput "fileOutput:$binFile,"

		    ## number of data points
		    set Ndata [ expr int($duration * $samplingF) ]

		    ## amount to skip for initial transient
		    set sliceStart [ expr int($TransientSkip * $samplingF) ]

		    ## length of psd
		    set psdlength [ expr $Ndata / 16 ]
		    
		    ## end time of segment
		    set etime [ expr $start_time + $duration + $TransientSkip ]

		    ## times in format start-end
		    set times $start_time-$etime

		    ## channel name alias
		    set chalias [ mkchannelnt $channel ]

		    ## numerical copy
		    set cst [ expr $start_time ]
		    set cet [ expr $etime ]
		    set ctimes $cst-$cet



		    ################################################################################
		    ## Calibration
		    ################################################################################

		    if { $start_time < 734400013 } {
			
			## S2

			## Frame type for CalFac
			set FTypeCALFAC CAL_FAC_V03
			append FTypeCALFAC "_$IFO2"


			## L1 reference times for calibration
			set l1calstart 731488397

			## H1 reference times for calibration
			set h1calstart 734073939

			## H2 reference times for calibration
			if { $start_time <= 731849042 } {
			    set h2calstart 734234126
			} else {
			    set h2calstart 734234127
			}

			## Times for reference data
			set h1caltimes ${h1calstart}-[expr $h1calstart + 64 ]
			set h2caltimes ${h2calstart}-[expr $h2calstart + 64 ]
			set l1caltimes ${l1calstart}-[expr $l1calstart + 64 ]

			## Get IFO independent time
			set caltimes $ifo2
			append caltimes "caltimes" 
			set caltimes [ set $caltimes ]
			
			## set channel names
			set h1oloop [ mkchannelnt $IFO2:CAL-OLOOP_FAC ]
			set h1cavfac [ mkchannelnt $IFO2:CAL-CAV_FAC ]
			set h1gain [ mkchannelnt $IFO2:CAL-CAV_GAIN ]
			set h1resp [ mkchannelnt $IFO2:CAL-RESPONSE ]


			if { $IFO2 == "H2" } {
			    ## For H2, hard-coded location of reference calibration
			    set frqueryREF "CAL_REF_V03_$IFO2 $IFO /ldas_outgoing/mirror/frames/S2/LHO/cal/H-CAL_REF_V03_H2-${h2calstart}-64.gwf $caltimes proc($IFO2:CAL-CAV_GAIN!0!7000.0001!) h1gain\nCAL_REF_V03_$IFO2 $IFO /ldas_outgoing/mirror/frames/S2/LHO/cal/H-CAL_REF_V03_H2-${h2calstart}-64.gwf $caltimes proc($IFO2:CAL-RESPONSE!0!7000.0001!) h1resp"
			} else {
			    ## For H1 or L1, used standard file location
			    set frqueryREF "CAL_REF_V03_$IFO2 $IFO $caltimes proc($IFO2:CAL-CAV_GAIN!0!7000.0001!) h1gain\nCAL_REF_V03_$IFO2 $IFO $caltimes proc($IFO2:CAL-RESPONSE!0!7000.0001!) h1resp "
			} 

			## Frame query for calibration frames
			set frqueryFAC "$FTypeCALFAC $IFO $ctimes proc($IFO2:CAL-OLOOP_FAC) h1oloop\n$FTypeCALFAC $IFO $ctimes proc($IFO2:CAL-CAV_FAC) h1cavfac"
			
			
			## Calibration aliases (LDAS only)
			set calaliases "                   h1oloop  = $h1oloop;
                   h1cavfac = $h1cavfac;
                   h1gain   = $h1gain;
                   h1resp   = $h1resp;"

			## Calibration algorithms
			set calalgo "
                        cavfac = float(h1cavfac);
                        cavfaccplx = complex(cavfac);

                        output(cavfaccplx,_,_,$IFO2:CAL-CAV_FAC,$IFO2 cavity factor \[COMPLEX8TimeSeries\]);

                        oloop = float(h1oloop);
                        oloopcplx = complex(oloop);

                        output(oloopcplx,_,_,$IFO2:CAL-OLOOP_FAC,$IFO2 open loop factor \[COMPLEX8TimeSeries\]);
                        output(h1gain,_,_,$IFO2:CAL-CAV_GAIN,$IFO2 reference cavity gain \[COMPLEX8FrequencySeries\]);
                        output(h1resp,_,_,$IFO2:CAL-RESPONSE,$IFO2 reference response \[COMPLEX8FrequencySeries\]);

                        h1cavfacf = float(h1cavfac);
                        h1cavfacf = complex(h1cavfacf);

                        h1oloopf = float(h1oloop);
                        h1oloopf = complex(h1oloopf);
"
		    }

		    
		    if { $start_time > 734400013 } {

			## S3

			set FTypeCALFAC CAL_FAC_V01
			append FTypeCALFAC "_$IFO2"

			# Magic times for calibration
			set l1caltimes "753424982-753425045"

			set h1caltimes "757806384-757806447"
			
			set h2caltimes "758175883-758175946"

			## get caltimes for the right IFO
			set caltimes $ifo2
			append caltimes "caltimes" 
			set caltimes [ set $caltimes ]

			## set ref frame type
			if { $IFO2 == "L1" } {
			    set CALREFTYPE "CAL_REF_V01P1_L1"
			} else {
			    set CALREFTYPE "CAL_REF_V01_$IFO2"
			}

			## define frame query
			set frqueryREF "$CALREFTYPE $IFO $caltimes proc($IFO2:CAL-CAV_GAIN!0!7000.0001!) h1gain\n$CALREFTYPE $IFO $caltimes proc($IFO2:CAL-RESPONSE!0!7000.0001!) h1resp"
		     

			## get alpha's/beta's
			set frqueryFAC "$FTypeCALFAC $IFO $ctimes proc($IFO2:CAL-OLOOP_FAC h1oloop\n$FTypeCALFAC $IFO $ctimes proc($IFO2:CAL-CAV_FAC) h1cavfac"


			## datacond algorithms for calibration
			set calalgo "
                        cavfac = float(h1cavfac);
                        cavfaccplx = complex(cavfac);

                        output(cavfaccplx,_,_,$IFO2:CAL-CAV_FAC,$IFO2 cavity factor \[COMPLEX8TimeSeries\]);

                        oloop = float(h1oloop);
                        oloopcplx = complex(oloop);

                        output(oloopcplx,_,_,$IFO2:CAL-OLOOP_FAC,$IFO2 open loop factor \[COMPLEX8TimeSeries\]);
                        output(h1gain,_,_,$IFO2:CAL-CAV_GAIN,$IFO2 reference cavity gain \[COMPLEX8FrequencySeries\]);
                        output(h1resp,_,_,$IFO2:CAL-RESPONSE,$IFO2 reference response \[COMPLEX8FrequencySeries\]);

                        h1cavfacf = float(h1cavfac);
                        h1cavfacf = complex(h1cavfacf);

                        h1oloopf = float(h1oloop);
                        h1oloopf = complex(h1oloopf);
"
		    }


		    ## if NoCalibration is defined:
		    if { [ info exists NoCalibration ] } {
			set frqueryREF ""
			set frqueryFAC ""
			set calaliases ""
			set calalgo ""
		    }


		    ## If we need to do MDC injections:
		    if { [ info exists MDCFrames ] } {
			## frame query
			set frqueryMDC "$MDCFrames HL $times adc(${IFO2}:GW) y"

			## algorithms to add to data
			set MDCalias "y = ${IFO2}\\:GW"
			set MDCalgo "x = add(x0,y);"
		    } else {
			## No MDC injections
			set frqueryMDC ""
			set MDCalgo "x = value(x0);"
			set MDCalias ""
		    }

		    ################################################################################
		    ## Job submission
		    ################################################################################

		    ## manager
		    if { $manager == "cit" } {
			set dataserver "ldas-gridmon.ligo.caltech.edu"
		    } else {
			error "Unknown manager $manager"
		    }

		    ## frames query and aliases 
		    set framequery "$FType $IFO $times adc($channel) x0\n$frqueryFAC\n$frqueryREF\n$frqueryMDC"

		    ## algorithms
		    set algorithms "$MDCalgo
                gwchn = double(x);
                $dcfilters
                gwchns = slice(gwchn,$sliceStart,$Ndata,1);
                gwchns = float(gwchns);
                output(gwchns,_,_,GW_STRAIN_DATA:primary,GW_STRAIN_DATA);
                $calalgo
                $dcwave
                spec = psd( gwchns, $psdlength );
                output(spec,_,_,GW_STRAIN_PSD,GW_STRAIN_PSD);
"

		    ## filterParams
		    set filterparams "-filterparams,$binoutput$burstoutput$channel,$Ndata,$injWave,$injAmp,$injAlpha,$injDelta,$injPsi,$injN,$injTimes,$ETG,$ETGParams"

		    ## response files
		    set responsefiles $responseFile

		    ## Check we have a valid proxy
		    if { [ catch { exec grid-proxy-info } cout ] } {
			puts $cout
		    } else {
			if { [ regexp {timeleft : ([0-9]+:[0-9]+:[0-9]+)} $cout junk timleft ] == 0} {
			    error "Invalid output from grid-proxy-info: $cout"
			    exit
			} else {

			    regexp {([0-9]+):([0-9]+):([0-9]+)} $timleft junk h m s 

			    set h [ string trimleft $h 0 ]
			    if { [ string length $h ] == 0 } {
				set h 0
			    }

			    set m [ string trimleft $m 0 ]
			    if { [ string length $m ] == 0 } {
				set m 0
			    }

			    set s [ string trimleft $s 0 ]
			    if { [ string length $s ] == 0 } {
				set s 0
			    }

			    set tup [ expr 3600*int($h)+60*int($m)+int($s) ]

			    if { $tup > 3600 } {
				puts "proxy is up ($tup s left)"
			    } else {
				puts "Less than 1 hour left: need to initialize proxy"
				puts "Please run grid-proxy-init"
			    }
			}
		    }

		    ## write framequery.txt
		    set fid [ open "framequery.txt" "w" ]
		    puts $fid $framequery
		    close $fid

		    ## write algorithms.txt
		    set fid [ open "algorithms.txt" "w" ]
		    puts $fid $algorithms
		    close $fid

		    ## write filterparams.txt
		    set fid [ open "filterparams.txt" "w" ]
		    puts $fid $filterparams
		    close $fid

		    ## write responsefiles.txt
		    set fid [ open "responsefiles.txt" "w" ]
		    puts $fid $responsefiles
		    close $fid

		    set jobOK 0

	   
		    if { [ info exists NoCondor ] } {
			## run outside of Condor (stand-alone)
			if { [ catch { exec burstdso $dataserver framequery.txt algorithms.txt filterparams.txt responsefiles.txt 2>> std.err } cout ] } {
			    puts $cout

			    if { $fileOutput == 1 } {
				if { [ catch { exec mv $binFile "$binFile.failed" } cout ] } {
				    puts $cout
				}
			    }

			    set eti [ expr $start_time + $duration ]
			    puts "Giving up for $start_time $eti"
			    incr start_time $duration

			    ## clean up
			    if { [ info exists injType ] } {
				if { $injType >= 2 } {
				    catch { exec rm hamp.ilwd.$start_time }
				    catch { exec rm alpha.ilwd.$start_time }
				    catch { exec rm delta.ilwd.$start_time }
				    catch { exec rm psi.ilwd.$start_time }
				    catch { exec rm injTime.ilwd.$start_time }

				    catch { exec rm hamp${IFO2}.ilwd.$start_time }
				    catch { exec rm alpha${IFO2}.ilwd.$start_time }
				    catch { exec rm delta${IFO2}.ilwd.$start_time }
				    catch { exec rm psi${IFO2}.ilwd.$start_time }
				    catch { exec rm injTime${IFO2}.ilwd.$start_time }
				}
			    }
			    
			} else {
			    set jobOK 1
			    set jobRetry 0
			}



			if { $jobOK } {
			    
			    if { $moveIt } {
				set zeFile ./results$channel$Job1(jobid).bin
				
				if { [ catch { exec mv $binFile $zeFile } cout ] } {
				    puts $cout
				}
			    } else {
				set zeFile $binFile
			    }
			    
			}
			
			if { $doSplit } {
			    if { [ catch { exec sblit $zeFile >& toc$IFO2.$Job1(jobid).$bId } cout ] } {
				puts $cout
			    } else {
				file delete $zeFile
			    }
			} else {
			    
			    if { $Zip } {
				
				if { [ catch { exec gzip $zeFile } cout ] } {
				    puts $cout
				}
				
			    }
			    
			}


		    } else {
			## run with Condor
			puts "Creating condor job..."

			## copy files to working dir
			file copy -force algorithms.txt framequery.txt filterparams.txt responsefiles.txt $PrebinFile

			## give them a unique name
			file rename -force $PrebinFile/algorithms.txt $PrebinFile/algorithms.$start_time.txt
			file rename -force $PrebinFile/framequery.txt $PrebinFile/framequery.$start_time.txt
			file rename -force $PrebinFile/filterparams.txt $PrebinFile/filterparams.$start_time.txt
			file rename -force $PrebinFile/responsefiles.txt $PrebinFile/responsefiles.$start_time.txt

			## tell condor to copy files to cluster nodes
			if { [ string length $CondorFiles ] > 0 } {
			    set CondorFiles "$CondorFiles,algorithms.$start_time.txt,framequery.$start_time.txt,filterparams.$start_time.txt,responsefiles.$start_time.txt"
			} else {
			    set CondorFiles "algorithms.$start_time.txt,framequery.$start_time.txt,filterparams.$start_time.txt,responsefiles.$start_time.txt"
			}

			## create condor submit file
			if { $cjid == 0 } {
			    exec echo "\#!/bin/bash\n\n" > tmp.condor
			    exec echo "export PATH=$env(PATH):\$PATH" >> tmp.condor
			    exec echo "export LD_LIBRARY_PATH=$env(LD_LIBRARY_PATH):\$LD_LIBRARY_PATH" >> tmp.condor
			    exec echo "export BURSTWRAPPER=$env(BURSTWRAPPER)" >> tmp.condor
			    exec echo "export X509_USER_PROXY=$env(X509_USER_PROXY)\n\n" >> tmp.condor

			    if { [ file exists tmp ] } {
				file delete tmp
			    }

			    exec echo "universe = vanilla\n\nexecutable = IDIR/EXEC\n\nerror = job.EXEC.err\nlog = jobs.log\n\n\ninitialdir = IDIR\n\n" > tmp

			} 

			## Only query file location once per lock segment
			if { $EndSegTime < $start_time } {

			    set EndSegTime $end_time

			    ## create framequery for full segment

			    ## Parse framequery.txt
			    set fid [ open $PrebinFile/framequery.$start_time.txt "r" ]
			    while { [ eof $fid ] == 0 } {
				gets $fid line

				## get times for GW data
				if { [ regexp "$FType $IFO (\[0-9\]{9})-(\[0-9\]{9}) .*" $line junk T0 T1 ] } {
				    break
				}
			    }
			    close $fid

			    ## substitute T0 with start_time in $PrebinFile/framequery.$start_time.txt, write 
			    ## result in sfquery.txt
			    if { [ catch { exec subst $PrebinFile/framequery.$start_time.txt sfquery.txt $T0 $start_time } cout ] } {
				puts $cout
				exit 1
			    }
    
			    ## replace T1 with end_time
			    if { [ catch { exec subst sfquery.txt $PrebinFile/sframequery.txt $T1 $end_time } cout ] } {
				puts $cout
				exit 1
			    }
    
			    ## run buildcache utility on sframequery.txt
			    if { [ catch { exec /bin/bash --login -c "echo $PrebinFile ; cd $PrebinFile ; buildcache $dataserver $PrebinFile/sframequery.txt" } cout ] } {
				puts $cout
				exit 1
			    }

			    ## rename cache file
			    file copy -force $PrebinFile/FrCacheFile $PrebinFile/FrCacheFile.$start_time
			    
			    set FrCFile FrCacheFile.$start_time

			    ## tell condor to send cache file to nodes
			    set CondorFiles "$CondorFiles,FrCacheFile.$start_time"

			    set SegTime $start_time
			} else {
			    if { $cjid == 0 } {

				set CondorFiles "$CondorFiles,FrCacheFile.$SegTime"
				
			    }
			}

			## add burstdso command to tmp.condor file
			exec echo "burstdso $FrCFile framequery.$start_time.txt algorithms.$start_time.txt filterparams.$start_time.txt responsefiles.$start_time.txt\n" >> tmp.condor

			incr cjid


			## Check that we have enough jobs to submit efficiently to Condor
			if { $cjid > $NCondorJobs } {
			    puts "*********************************************************************"
			    puts "*********************************************************************"
			    puts "Sending jobs to cluster"
			    puts "*********************************************************************"
			    puts "*********************************************************************"

			    ## rename submit file
			    file rename -force tmp.condor $PrebinFile/tmp.condor.$SegTime

			    ## which files to transfer to nodes
			    exec echo "transfer_input_files = $CondorFiles" >> tmp


			    exec echo "when_to_transfer_output = ON_EXIT\n\ngetenv = True\n\nQueue\n\n" >> tmp

			    ## replace EXEC with tmp.condor.$SegTime
			    if { [ catch { exec subst tmp tmp2 EXEC tmp.condor.$SegTime } cout ] } {
				puts $cout
				exit 1
			    }

			    ## escape slashes
			    regsub -all "/" $ZeFile "\\/" rZeFile

			    ## replace IDIR with $rZeFile
			    if { [ catch { exec subst tmp2 tmp3 IDIR $rZeFile } cout ] } {
				puts $cout
				exit 1
			    }

    
			    ## call condor_submit
			    if { [ catch { exec /bin/bash --login -c "condor_submit tmp3" } cout ] } {
				puts $cout
				exit 1
			    }


			    ## clean up
			    file delete tmp.condor

			    set cjid 0
			    set CondorFiles ""

			}



		    }


		    ## end loop over channels:
		} 
    
		## got to next segment
		incr start_time $duration

		## end if do1:
	    }

	    ## end loop over times within one segment:
	}

	## end if gotIFO, etc.
    }

    ## end loop over segments:
}


## Run any left over jobs
if { [ info exists NoCondor ] == 0 && $cjid > 0 } {
    puts "*********************************************************************"
    puts "*********************************************************************"
    puts "Sending jobs to cluster"
    puts "*********************************************************************"
    puts "*********************************************************************"

    file rename -force tmp.condor $PrebinFile/tmp.condor.$SegTime
    
    exec echo "transfer_input_files = $CondorFiles" >> tmp


    exec echo "when_to_transfer_output = ON_EXIT\n\ngetenv = True\n\nQueue\n\n" >> tmp


    if { [ catch { exec subst tmp tmp2 EXEC tmp.condor.$SegTime } cout ] } {
	puts $cout
	exit 1
    }

    regsub -all "/" $ZeFile "\\/" rZeFile
    if { [ catch { exec subst tmp2 tmp3 IDIR $rZeFile } cout ] } {
	puts $cout
	exit 1
    }

    

    if { [ catch { exec /bin/bash --login -c "condor_submit tmp3" } cout ] } {
	puts $cout
	exit 1
    }


    file rename -force tmp3 tmp3.$SegTime

    set cjid 0
    set CondorFiles ""
}


puts "Normal exit"

## the end.
