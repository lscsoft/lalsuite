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

proc mkchannel { chan stime } {
    regsub -all -- {:} $chan {\\:} chalias
    set chalias ${chalias}::AdcData:$stime:0:Frame
}

proc mkchannelnt { chan } {
    regsub -all -- {:} $chan {\:} chalias
    return $chalias
}

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

################################
set injAmpTmp $injAmp
set LOCALFILE 0
set FastCacheDone 0

set EndSegTime 0
set cjid 0
set CondorFiles ""

if { [ regexp {LOCAL:(.*)} $PrebinFile junk ZeFile ] } {

    set NNODES 210

    set PrebinFile $ZeFile

    if { [ file exists $PrebinFile ] == 0 } {
	file mkdir $PrebinFile
    }

    if { [ regexp {/usr1/(.*)} $ZeFile junk PrebinFileSUFFIX ] == 0 } {
	error "PrebinFile doesn't start with /usr1/"
    }

    for { set ni 1 } { $ni < $NNODES } { incr ni 1 } {

	if { [ file exists "/data/node${ni}" ] } {
	if { [ file exists "/data/node${ni}/$PrebinFileSUFFIX" ] == 0 } {
	    puts "Creating /data/node${ni}/$PrebinFileSUFFIX"
	    file mkdir "/data/node${ni}/$PrebinFileSUFFIX"

	    if { [ catch { exec chmod 777 "/data/node${ni}/$PrebinFileSUFFIX" } err ] } {
		puts "WARNING: can't change mode of /data/node${ni}/$PrebinFileSUFFIX to 777; LDAS might not be able to write there!!"
	    }
	} else {
	    puts "WARNING: /data/node${ni}/$PrebinFileSUFFIX already exists!"
	}
	}
    }

} else {

    if { [ file exists $PrebinFile ] == 0 } {
	puts "Creating $PrebinFile"
	file mkdir $PrebinFile
    } else {
	puts "WARNING: $PrebinFile already exists!"
    }

    if { [ catch { exec chmod 777 $PrebinFile } err ] } {
	puts "WARNING: can't change mode of $PrebinFile to 777; LDAS might not be able to write there!!"
    }

}

set samplingF [ expr int($samplingF) ]

set channels [ split $channel "," ]

set IFOs [ list ]
foreach channel $channels {
    lappend IFOs [ string range $channel 0 1 ]
}

set fileOutput 0
set doSplit 0

switch -exact -- $OutputType {
    0 { set binputput "" }
    1 { set binoutput "binOutput," ; set doSplit 1 }
    2 { set fileOutput 1 }
    default { error "Invalue OutputType: $OutputType" }
}

if { $BurstOutput == 0 } {
    set burstoutput ""
} else {
    set burstoutput "noLALBurstOutput,"
}

if { [ llength $argv ] == 2 } {
    set doblock 1
    set blockId [ lindex $argv 1 ]
} else {
    set doblock 0
    set blockId 0
}

eval set dataType $dataType

switch -exact -- $dataType {
    0 { 
	set playgndURL "http://www.ligo.caltech.edu/~jsylvest/S2/coin/Playgnd.txt" 
	set SegmentList [ dataDownload $playgndURL ]
	set SegmentList [ split $SegmentList "\n" ]
    }
    1 { error "Production data mode not yet supported" }
    default { 
#	set SegmentFile $dataType
#	if [catch {open $SegmentFile r+} jID] {
#	    puts "Cannot open $SegmentFile for reading: $jID"
#	    exit 1
#	}    
#	set SegmentList ""
#	while {[gets $jID line] >= 0} {
#	    append SegmentList $line \n
#	}
#	close $jID
#	set SegmentList [ split $SegmentList "\n" ]

	set NSL 0
	set SegmentFileList [ split $dataType "," ]
	foreach file $SegmentFileList {
	    set SegmentListTmp [SegRead $file]
	    if { $NSL < 1 } {
		set SegmentList $SegmentListTmp
	    } else {
		set SegmentList [SegIntersection $SegmentList $SegmentListTmp]
	    }
	    
	    incr NSL 1
	}
    }
}

set jobRetry 0
set bId 0 
set do1 1


foreach l1 $SegmentList {

    if { $dataType <= 1 } {
	if { [ regexp {([0-9\.]+)\s+([0-9\.]+)\s+([01])\s+([01])\s+([01]).*} $l1 junk start_time end_time H1 H2 L1 ] == 0 } {
	    continue
	}
    } else {
	set start_time [lindex $l1 0]
	set end_time   [lindex $l1 1]  
	incr end_time  1

	foreach I $IFOs {
#	    if { [ regexp $I $dataType ] } {
		set $I 1
#	    } else {
#		set $I 0
#	    }
	}
    }

    set start_time [ expr int($start_time) ]

    incr bId 1

    set gotIFO 0
    foreach I $IFOs {
	set gotIFO [ set $I ]
	if { $gotIFO == 0 } {
	    break
	}
    }

  if { $gotIFO && ( $doblock == 0 || $bId >= $blockId ) } {

  while { [ expr $start_time + $duration ] <= $end_time } {

    set seed_return [ expr srand($start_time) ]


    if { $do1 == $do1inN } {
	set do1 0
    }
    incr do1 1

    if { $do1 == 1 } {

	set channelNumber 0

	foreach channel $channels {

	    incr channelNumber 1

	    set IFO [ string index $channel 0 ]
	    set IFO2 [ string range $channel 0 1 ]
	    set ifo2 [ string tolower $IFO2 ]

    set PbId [SegID S2 $IFO2 $start_time]

    puts "Channel $channel: processing $start_time for $duration seconds (block $bId - Peter's block Id $PbId)"


################################
if { [ string length $prefilters ] == 0 } {
    set dcfilters ""
    set responseFile ""
} else {
    set flist [ split $prefilters "," ]
    if { [ llength $flist ] == 0 } {
	error "Invalid filter list: $prefilters"
    }

    set responseFile ""
    set dcfilters ""
    set fcount 1
    set got1 0
    set fa [ lindex $flist 0 ]
    append fa "_a.ilwd"
    if { [ file exists $fa ] } {
	append responseFile "$fa push fa$fcount\n"
	incr got1 1
    } else {
	append responseFile "${filtroot}/$fa push fa$fcount\n"
	incr got1 1
    }
    set fb [ lindex $flist 0 ]
    append fb "_b.ilwd"
    if { [ file exists $fb ] } {
	append responseFile "$fb push fb$fcount\n"
	incr got1 1
    } else {
	append responseFile "${filtroot}/$fb push fb$fcount\n"
	incr got1 1
    }
    if { $got1 != 2 } {
	error "Can't file $fa or $fb"
    }
    append dcfilters "\t\tgwchn = linfilt(fb$fcount,fa$fcount,gwchn);\n"

    incr fcount 1

    if { [ llength $flist ] > 1 } {

	foreach f [ lrange $flist 1 end ] {

	    set got1 0
	    set fa $f
	    append fa "_a.ilwd"
	    if { [ file exists $fa ] } {
		append responseFile "$fa push fa$fcount\n"
		incr got1 1
	    } else {
		append responseFile "${filtroot}/$fa push fa$fcount\n"
		incr got1 1
	    }
	    set fb $f
	    append fb "_b.ilwd"
	    if { [ file exists $fb ] } {
		append responseFile "$fb push fb$fcount\n"
		incr got1 1
	    } else {
		append responseFile "${filtroot}/$fb push fb$fcount\n"
		incr got1 1
	    }
	    if { $got1 != 2 } {
		error "Can't file $fa or $fb"
	    }
	    append dcfilters "\t\tgwchn = linfilt(fb$fcount,fa$fcount,gwchn);\n"

	    incr fcount 1

	}

    }
}


if { [ string length $prefiltersB ] > 0 } {

    set flist [ split $prefiltersB "," ]
    set got1 0
    set fb [ lindex $flist 0 ]

    eval set fb $fb 

    append fb ".ilwd"
    if { [ file exists $fb ] } {
	append responseFile "$fb push fb$fcount\n"
	incr got1 1
    } else {
	append responseFile "${filtroot}/$fb push fb$fcount\n"
	incr got1 1
    }
    append dcfilters "\t\tgwchn = linfilt(fb$fcount,gwchn);\n"
    incr fcount 1
    if { [ llength $flist ] > 1 } {
	foreach f [ lrange $flist 1 end ] {
	    set got1 0
	    set fb $f

	    eval set fb $fb 
	    
	    append fb ".ilwd"
	    if { [ file exists $fb ] } {
		append responseFile "$fb push fb$fcount\n"
		incr got1 1
	    } else {
		append responseFile "${filtroot}/$fb push fb$fcount\n"
		incr got1 1
	    }
	    append dcfilters "\t\tgwchn = linfilt(fb$fcount,gwchn);\n"
	    incr fcount 1
	}
    } 
}

regsub -all "gwchn" $dcfilters "hpf" hpfilters
regsub -all "gwchn" $dcfilters "hcf" hcfilters


if { [ string length $waveforms ] == 0 } {
    set injWave ""
    set dcwave ""
} else {
    set wlist [ split $waveforms "," ]
    if { [ llength $wlist ] == 0 } {
	error "Invalid waveform list: $waveforms"
    }

    set wcount 1
    set got0 0
    set injWave ""
    set dcwave ""
    
    set got1 0
    set got2 0
    set ff [ lindex $wlist 0 ]
    set fplus [ lindex $wlist 0 ]
    append fplus "_p.ilwd"
    if { [ file exists $fplus ] } {
	append responseFile "$fplus push ihp$wcount\n"
	set got1 1
    }
    set fcross [ lindex $wlist 0 ]
    append fcross "_c.ilwd"
    if { [ file exists $fcross ] } {
	append responseFile "$fcross push ihc$wcount\n"
	set got2 1
    }
    if { $got1 == 0 && $got2 == 0 } {
	error "Can't file $fplus or $fcross"
    }
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
	append dcwave "\t\tihpt = tseries(ihp$wcount,$samplingF.0,$start_time);\n"
	append dcwave "\t\thpf = respfilt(ihpt,h1resp,h1gain,h1cavfacf,h1oloopf);\n"
	append dcwave $hpfilters
	append dcwave "\t\thpfr = float(hpf);\n"
	set st1 "\t\toutput(hpfr,_,_,$ff"
	append dcwave [ append st1 "_p,$ff" "_p);\n" ]
    } else {
	set st1 "\t\toutput(zero,_,_,$ff"
	append dcwave [ append st1 "_p,$ff" "_p);\n" ]
    }
    if { $got2 != 0 } {
	append dcwave "\t\tihct = tseries(ihc$wcount,$samplingF.0,$start_time);\n"
	append dcwave "\t\thcf = respfilt(ihct,h1resp,h1gain,h1cavfacf,h1oloopf);\n"
	append dcwave $hcfilters
	append dcwave "\t\thcfr = float(hcf);\n"
	set st1 "\t\toutput(hcfr,_,_,$ff"
	append dcwave [ append st1 "_c,$ff" "_c);\n" ]
    } else {
	set st1 "\t\toutput(zero,_,_,$ff"
	append dcwave [ append st1 "_c,$ff" "_c);\n" ]
    }
    incr wcount 1
    

    set w _[ lindex $wlist 0 ]
    set injWave [ append w "_" ]

    if { [ llength $wlist ] > 1 } {
	set injWave "($injWave"
	foreach w [ lrange $wlist 1 end ] {

	    append injWave ",_" $w "_"



	    set got1 0
	    set got2 0
	    set fplus $w
	    append fplus "_p.ilwd"
	    if { [ file exists $fplus ] } {
		append responseFile "$fplus push ihp$wcount\n"
		set got1 1
	    }
	    set fcross $w
	    append fcross "_c.ilwd"
	    if { [ file exists $fcross ] } {
		append responseFile "$fcross push ihc$wcount\n"
		set got2 1
	    }
	    if { $got1 == 0 && $got2 == 0 } {
		error "Can't file $fplus or $fcross"
	    }
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
		append dcwave "\t\tihpt = tseries(ihp$wcount,$samplingF.0,$start_time);\n"
		append dcwave "\t\thpf = respfilt(ihpt,h1resp,h1gain,h1cavfacf,h1oloopf);\n"
		append dcwave $hpfilters
		append dcwave "\t\thpfr = float(hpf);\n"
		set st1 "\t\toutput(hpfr,_,_,$w"
		append dcwave [ append st1 "_p,$w" "_p);\n" ]
	    } else {
		set st1 "\t\toutput(zero,_,_,$w"
		append dcwave [ append st1 "_p,$w" "_p);\n" ]
	    }
	    if { $got2 != 0 } {
		append dcwave "\t\tihct = tseries(ihc$wcount,$samplingF.0,$start_time);\n"
		append dcwave "\t\thcf = respfilt(ihct,h1resp,h1gain,h1cavfacf,h1oloopf);\n"
		append dcwave $hcfilters
		append dcwave "\t\thcfr = float(hcf);\n"
		set st1 "\t\toutput(hcfr,_,_,$w"
		append dcwave [ append st1 "_c,$w" "_c);\n" ]
	    } else {
		set st1 "\t\toutput(zero,_,_,$w"
		append dcwave [ append st1 "_c,$w" "_c);\n" ]
	    }
	    incr wcount 1

	}
	append injWave ")"
    }
}


if { [ info exists injType ] == 1 } {
    switch -exact -- $injType {
	1 { 

	    if { $channelNumber == 1 } {
		if { [ catch { exec ShellInjection1 $injAmpTmp $injNtot $injN $duration } cout ] } {
		    puts $cout
		    exit -1
		}
	    }

	    set injAmp "__hamp__"
	    set injAlpha "__alpha__"
	    set injDelta "__delta__"
	    set injPsi "__psi__"
	    set injTimes "__injTime__"

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

	    if { $channelNumber == 1 } {
		if { [ catch { exec ShellInjection $injAmpTmp $injNtot $injN $duration } cout ] } {
		    puts $cout
		    exit -1
		}
	    }

	    set injAmp "__hamp__"
	    set injAlpha "__alpha__"
	    set injDelta "__delta__"
	    set injPsi "__psi__"
	    set injTimes "__injTime__"

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

	    set skip 0
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

	    set injAmp "__hamp${IFO2}__"
	    set injAlpha "__alpha${IFO2}__"
	    set injDelta "__delta${IFO2}__"
	    set injPsi "__psi${IFO2}__"
	    set injTimes "__injTime${IFO2}__"

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


################################

if { [ info exists ETGParameters ] == 0 } { 
    if { $ETG == "TFCLUSTERS" } {
	set TFCchn [ set channel ]
	set ETGParams "$TFCchn,$TFCThr,$TFCWin,$TFCThrMethod,$TFCdata,$TFCTRez,$TFCFmin,$TFCFmax,$TFCalpha,$TFCsigma,$TFCdelta"
    }
} else {
    eval set ETGParameters $ETGParameters

    if { [ regexp {__(.+)__} $ETGParameters junk file ] } {
	append responseFile "${file}.ilwd pass\n"
    }
    set ETGParams $ETGParameters
}




#    set now [ clock seconds ]

    set binFile "$PrebinFile/job$IFO2.$userName.$start_time.$bId.bin"
    set binoutput "fileOutput:$binFile,"

    set Ndata [ expr int($duration * $samplingF) ]
    set sliceStart [ expr int($TransientSkip * $samplingF) ]

    set psdlength [ expr $Ndata / 16 ]
    set etime [ expr $start_time + $duration + $TransientSkip ]
    set times $start_time-$etime
    #set chalias [ mkchannel $channel $start_time ]
    set chalias [ mkchannelnt $channel ]

    set cst [ expr $start_time ]
    set cet [ expr $etime ]
    set ctimes $cst-$cet

    #set FTypeCALFAC "CAL_FAC"
    #set FTypeCALFAC $FType
#    set FTypeCALFAC SenseMonitor_$IFO2
#    append FTypeCALFAC "_M" 

    set FTypeCALFAC CAL_FAC_V03
    append FTypeCALFAC "_$IFO2"

    # Magic times for calibration
    set l1calstart 731488397
    set h1calstart 734073939

    if { $start_time <= 731849042 } {
	set h2calstart 734234126
    } else {
	set h2calstart 734234127
    }

    set h1caltimes ${h1calstart}-[expr $h1calstart + 64 ]
    set h2caltimes ${h2calstart}-[expr $h2calstart + 64 ]
    set l1caltimes ${l1calstart}-[expr $l1calstart + 64 ]

    set caltimes $ifo2
    append caltimes "caltimes" 
    set caltimes [ set $caltimes ]

    set h1oloop [ mkchannelnt $IFO2:CAL-OLOOP_FAC ]
    set h1cavfac [ mkchannelnt $IFO2:CAL-CAV_FAC ]
    set h1gain [ mkchannelnt $IFO2:CAL-CAV_GAIN ]
    set h1resp [ mkchannelnt $IFO2:CAL-RESPONSE ]


    if { $IFO2 == "H2" } {
	set frqueryREF "CAL_REF_V03_$IFO2 $IFO /ldas_outgoing/mirror/frames/S2/LHO/cal/H-CAL_REF_V03_H2-${h2calstart}-64.gwf $caltimes proc($IFO2:CAL-CAV_GAIN!0!7000.0001!) h1gain\nCAL_REF_V03_$IFO2 $IFO /ldas_outgoing/mirror/frames/S2/LHO/cal/H-CAL_REF_V03_H2-${h2calstart}-64.gwf $caltimes proc($IFO2:CAL-RESPONSE!0!7000.0001!) h1resp"
#	set frqueryREF "CAL_REF_V03_$IFO2 $IFO $caltimes proc($IFO2:CAL-CAV_GAIN!0!7000.0001!) h1gain\nCAL_REF_V03_$IFO2 $IFO $caltimes proc($IFO2:CAL-RESPONSE!0!7000.0001!) h1resp"
    } else {
	set frqueryREF "CAL_REF_V03_$IFO2 $IFO $caltimes proc($IFO2:CAL-CAV_GAIN!0!7000.0001!) h1gain\nCAL_REF_V03_$IFO2 $IFO $caltimes proc($IFO2:CAL-RESPONSE!0!7000.0001!) h1resp "
    } 

    set frqueryFAC "$FTypeCALFAC $IFO $ctimes proc($IFO2:CAL-OLOOP_FAC) h1oloop\n$FTypeCALFAC $IFO $ctimes proc($IFO2:CAL-CAV_FAC) h1cavfac"


    set calaliases "                   h1oloop  = $h1oloop;
                   h1cavfac = $h1cavfac;
                   h1gain   = $h1gain;
                   h1resp   = $h1resp;"

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

if { [ info exists $useResponse ] == 0 } {
    set useResponse ""
}

    if { [ info exists NoCalibration ] } {
	set frqueryREF ""
	set frqueryFAC ""
	set calaliases ""
	set calalgo ""
    }

if { [ info exists MDCFrames ] } {
    set frqueryMDC "$MDCFrames HL $times adc(${IFO2}:GW) y"

    set MDCalias "y = ${IFO2}\\:GW"
    set MDCalgo "x = add(x0,y);"

} else {
    set frqueryMDC ""
    set MDCalgo "x = value(x0);"
    set MDCalias ""
}

#######################################################################

# manager
if { $manager == "cit" } {
    set dataserver "ldas-gridmon.ligo.caltech.edu"
}

# frames query and aliases
set framequery "$frqueryFAC\n$frqueryREF\n$frqueryMDC"

# algorithm
set NNData [ expr int(($etime-$start_time)*16384) ]
set Seed [ clock clicks ]


set algorithms "
                seed = value($Seed);
                mone = sub(0,1);
                seed = mul(seed,mone);
                x1 = gasdev($NNData,seed);
                x0 = tseries(x1,16384.0,$start_time,0);
                $MDCalgo
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

# filterParams
set filterparams "-filterparams,$binoutput$burstoutput$channel,$Ndata,$injWave,$injAmp,$injAlpha,$injDelta,$injPsi,$injN,$injTimes,$ETG,$ETGParams"

# response files
set responsefiles $responseFile

#######################################################

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

set fid [ open "framequery.txt" "w" ]
puts $fid $framequery
close $fid

set fid [ open "algorithms.txt" "w" ]
puts $fid $algorithms
close $fid

set fid [ open "filterparams.txt" "w" ]
puts $fid $filterparams
close $fid

set fid [ open "responsefiles.txt" "w" ]
puts $fid $responsefiles
close $fid

set jobOK 0

if { [ info exists NoCondor ] } {
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


#######################################################

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
# Condor
puts "Creating condor job..."

file copy -force algorithms.txt framequery.txt filterparams.txt responsefiles.txt $PrebinFile

file rename -force $PrebinFile/algorithms.txt $PrebinFile/algorithms.$start_time.txt
file rename -force $PrebinFile/framequery.txt $PrebinFile/framequery.$start_time.txt
file rename -force $PrebinFile/filterparams.txt $PrebinFile/filterparams.$start_time.txt
file rename -force $PrebinFile/responsefiles.txt $PrebinFile/responsefiles.$start_time.txt

if { [ string length $CondorFiles ] > 0 } {
    set CondorFiles "$CondorFiles,algorithms.$start_time.txt,framequery.$start_time.txt,filterparams.$start_time.txt,responsefiles.$start_time.txt"
} else {
    set CondorFiles "algorithms.$start_time.txt,framequery.$start_time.txt,filterparams.$start_time.txt,responsefiles.$start_time.txt"
}

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

if { $EndSegTime < $start_time } {

    set EndSegTime $end_time

    # create framequery for full segment
    set fid [ open $PrebinFile/framequery.$start_time.txt "r" ]
    set gotFD 0
    while { [ eof $fid ] == 0 } {
	gets $fid line

	if { [ regexp "$FType $IFO (\[0-9\]{9})-(\[0-9\]{9}) .*" $line junk T0 T1 ] } {
	    set gotFD 1
	    break
	}
    }
    close $fid

    if { $gotFD == 0 } {
	set T0 $start_time
	set T1 $etime
    }

    if { [ catch { exec subst $PrebinFile/framequery.$start_time.txt sfquery.txt $T0 $start_time } cout ] } {
	puts $cout
	exit 1
    }
    
    if { [ catch { exec subst sfquery.txt $PrebinFile/sframequery.txt $T1 $end_time } cout ] } {
	puts $cout
	exit 1
    }
    

    if { [ catch { exec /bin/bash --login -c "echo $PrebinFile ; cd $PrebinFile ; buildcache $dataserver $PrebinFile/sframequery.txt" } cout ] } {
	puts $cout
	exit 1
    }

    file copy -force $PrebinFile/FrCacheFile $PrebinFile/FrCacheFile.$start_time

    set FrCFile FrCacheFile.$start_time

    set CondorFiles "$CondorFiles,FrCacheFile.$start_time"

    set SegTime $start_time
} else {
    if { $cjid == 0 } {

	set CondorFiles "$CondorFiles,FrCacheFile.$SegTime"
	
    }
}


exec echo "burstdso $FrCFile framequery.$start_time.txt algorithms.$start_time.txt filterparams.$start_time.txt responsefiles.$start_time.txt\n" >> tmp.condor

incr cjid



if { $cjid > $NCondorJobs } {
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





}


#end loop over channels:
} 
    
incr start_time $duration

# end if do1:
}

#end loop over times within one segment:
}

# end if gotIFO, etc.
}

# end loop over segments:
}


# Run any left over jobs
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
