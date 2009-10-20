#!/usr/bin/env tclsh

#MainSubmit.tcl - Submits the new dag to be started in the scheduler universe

foreach {var value} {
	help 0
	outlier_index unknown
	work_dir unknown
	} {
	global $var
	set $var [subst $value]
	}
	

proc usage {} {
puts {
Usage: MainSubmit.tcl [var value]

  help             display this message
  outlier_index    index of the outlier to be searched
  work_dir         working directory to run from
  
}
}

foreach {var value} $argv {
	global $var
	if {$var == "help"} {
		puts {
Usage: MainSubmit.tcl [var value]

  help             display this message
  outlier_index    index of the outlier to be searched
  work_dir         working directory to run from
  
			}
		}
	set $var $value
	}

cd $work_dir

source "configFile"

exec "condor_submit_dag $work_dir/idx$outlier_index/run_${iteration}.dag "

