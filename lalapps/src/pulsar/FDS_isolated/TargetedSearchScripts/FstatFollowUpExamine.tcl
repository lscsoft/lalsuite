#!/usr/bin/env tclsh


#ExamineFollowUpResults.py - Examines the results from a FStatFollowUp.py run, finds the largest events, makes a few plots, and calls FStatFollowUp.py if another iteration is still possible. 


set template_count_file "$outlier_dir/TemplateCount_${iteration}.txt"
exec $fstat_binary --Alpha [expr $ra - $ra_band/2.0 ] --Delta [expr $dec - $dec_band/2.0 ] --AlphaBand $ra_band  --DeltaBand $dec_band --Freq [expr $f0 - $f0_band/2.0 ] --FreqBand [expr $f0_band] --dFreq $freq_resolution --f1dot [expr $f1- $f1_band/2.0] --f1dotBand $f1_band --df1dot $dfreq_resolution --DataFiles ${sft_dir}/* --ephemDir $ephem_dir --ephemYear $ephem_year --refTime $param_time --gridType $grid_type --metricType $metric_type --metricMismatch $mismatch --countTemplates >> $template_count_file 2> /dev/null

set TEMPLATES_FILE [file open $template_count_file "r"]
while { ![eof $TEMPLATES_FILE] } {	
	gets $TEMPLATES_FILE line
	if {[lindex [split line ":"] 0] == "%% Number of templates"} {
		set total_templates [lindex [split line ":"] 1]
		}
	}

set sky_only_templates [expr $total_templates / ($f1_band/$dfreq_resolution)]

set total_templates [expr $total_templates * $f0_band/$freq_resolution] 

set angular_resolution = [expr sqrt($ra_band*$dec_band/$sky_only_templates)]



foreach {var value} {
	outlier_index unknown
	config_file unknown
	iteration unknown
	angular_resolution unknown
	total_templates unknown
	work_dir unknown
	} {
	global $var
	set $var [subst $value]
	}

foreach {var value} $argv {
	global $var
	set $var $value
	}
	
source config_file

###################################
#Checks to see if all the results files exist
#If some results files are missing, throw an error and exit

###################################
#Find the loudest

loudest2F = 0
loudestFreq = 0
loudestF1dot = 0
loudestF2dot = 0
loudestF3dot = 0
loudestRA = 0
loudestDEC = 0

array2F = array([])
arrayFreq = array([])
arrayRA = array([])
arrayDEC = array([])
arrayF1dot = array([])
arrayF2dot = array([])
arrayF3dot = array([])

recordLoudestResultFile = ''.join([Vars['current_directory'] + '/loudestResult_',str(Vars['iteration']),'.txt'])
record_loudest_results_file = open(recordLoudestResultFile,'w')

for job in range(0,Vars['jobs']):
    loudestFile = open(''.join([Vars['output_run_directory'],Vars['loudest_base_prefix'],str(job)]),'r')
    for line in loudestFile:
        if line.split('=')[0].strip() == 'twoF':
            temp2F = line.split('=')[1].strip().strip(';')
            
        elif line.split('=')[0].strip() == 'Alpha':
            tempRA = line.split('=')[1].strip().strip(';')
            
        elif line.split('=')[0].strip() == 'Delta':
            tempDEC = line.split('=')[1].strip().strip(';')
            
        elif line.split('=')[0].strip() == 'Freq':
            tempFreq = line.split('=')[1].strip().strip(';')

        elif line.split('=')[0].strip() == 'f1dot':
            tempF1dot = line.split('=')[1].strip().strip(';')

        elif line.split('=')[0].strip() == 'f2dot':
            tempF2dot = line.split('=')[1].strip().strip(';')

        elif line.split('=')[0].strip() == 'f3dot':
            tempF3dot = line.split('=')[1].strip().strip(';')

    loudestFile.close()

    record_loudest_results_file.write(' '.join([str(tempFreq),str(tempF1dot),str(tempF2dot),str(tempF3dot),str(tempRA),str(tempDEC),str(temp2F),"\n"]))

    if float(temp2F) > float(loudest2F):
        loudest2F = temp2F
        loudestRA = tempRA
        loudestDEC = tempDEC
        loudestFreq = tempFreq
        loudestF1dot = tempF1dot
        loudestF2dot = tempF2dot
        loudestF3dot = tempF3dot

#Finds all 2F values greater than or equal to a percentage (i.e. 90%) of the loudest 2F value, and
#takes a weighted average to find the most likely source position in frequency, spindown and sky position
#parameter space

average_cutoff = float(loudest2F) * float(Vars['percent_largest_to_average'])



for job in range(0,Vars['jobs']):
  fullResultsFile = open(''.join([Vars['output_run_directory'],Vars['results_base_prefix'],str(job)]),'r')
  for line in fullResultsFile:
    if not (line[0] == '%'):
      #Tests to see if 2F greater than required cutoff, then records data if true
      if float(line.split(' ')[6].strip()) >= average_cutoff:
        array2F= append(array2F,float(line.split(' ')[6].strip()))
        arrayFreq = append(arrayFreq, float(line.split(' ')[0].strip()))
        arrayRA = append(arrayRA,float(line.split(' ')[1].strip()))
        arrayDEC = append(arrayDEC,float(line.split(' ')[2].strip()))
        arrayF1dot = append(arrayF1dot,float(line.split(' ')[3].strip()))
        arrayF2dot = append(arrayF2dot,float(line.split(' ')[4].strip()))
        arrayF3dot = append(arrayF3dot,float(line.split(' ')[5].strip()))

#Calculate the weighted mean of all the parameters (weighted by 2F value)
total2Fsum = sum(array2F)
meanFreq = sum(array2F*arrayFreq)/total2Fsum
meanRA = sum(array2F*arrayRA)/total2Fsum
meanDEC = sum(array2F*arrayDEC)/total2Fsum
meanF1dot = sum(array2F*arrayF1dot)/total2Fsum
meanF2dot = sum(array2F*arrayF2dot)/total2Fsum
meanF3dot = sum(array2F*arrayF3dot)/total2Fsum

#Use the calculated average parameters to generate the parameter space for the next longer coherent step

Vars['Coherence_time'] = float(Vars['coherence_time'])
F_band = Vars['parameter_space_multiplier']*1.0/Vars['Coherence_time']
F_dot_band = Vars['parameter_space_multiplier']*1.0/Vars['Coherence_time']**2

RA_band = Vars['parameter_space_multiplier']*Vars['angular_resolution']
DEC_band = Vars['parameter_space_multiplier']*Vars['angular_resolution']

Vars['new_coherence_time'] = Vars['Coherence_time'] * 4.0


F_c = float(loudest2F)/2.0
Prob_below = (1 - ((1 + F_c)/2.0)*math.exp(-F_c))
FA = 1 - Prob_below**Vars['total_templates']




if os.path.exists(Vars['false_alarm_file']):
  false_alarm_file_exists = True
else:
  false_alarm_file_exists = False
  
FA_file = open(Vars['false_alarm_file'],'a')
if not false_alarm_file_exists:
  FA_file.write('Param_time Coherence_Time 2F RA DEC Freq F1dot F2dot F3dot Total-Templates FA Name Run\n')
FA_file.write(' '.join([str(Vars['param_time']),str(Vars['coherence_time']),loudest2F,loudestRA,loudestDEC,loudestFreq,loudestF1dot,loudestF2dot,loudestF3dot,str(Vars['total_templates']),str(FA),str(Vars['base_name']),str(Vars['iteration']),'\n']))
FA_file.close()

time.sleep(1)

#Cat result files for easy printing
os.system('cat ' + Vars['output_run_directory'].rstrip('/') + '/' + Vars['results_base_prefix'] + '* > ' + Vars['events_directory'].rstrip('/') + '/plot_data_' + Vars['run_name'])

#Delete the previous run to save space
if not Vars['messy']:
  os.system(''.join(['rm -rf ',Vars['base_name'],'_',str(Vars['iteration']),'_run']))

#Delete condor log files to save space
if not Vars['messy']:
  os.system(''.join(['rm node_',Vars['base_name'],'*']))


if FA < Vars['false_alarm_cutoff']:
  interesting_event_file = open(Vars['event_result_file'],'w')
  interesting_event_file.write('Param_time Coherence_Time 2F RA DEC Freq F1dot F2dot F3dot Total-Templates FA\n')
  interesting_event_file.write(' '.join([str(Vars['param_time']),str(Vars['coherence_time']),loudest2F,loudestRA,loudestDEC,loudestFreq,loudestF1dot,loudestF2dot,loudestF3dot,str(Vars['total_templates']),str(FA),'\n']))
  interesting_event_file.close()

  
  if float(Vars['Coherence_time']) >= (Vars['end_time'] - Vars['start_time']):
    print "Last search used all available data - no further iterations"
    sys.exit(0)


  #Command string for the next iteration to run
  command_string = ' '.join([''.join([Vars['python_dir'].rstrip('/'),'/FstatFollowUp.py']),
                             '-o',Vars['base_name'],
                             '-C',Vars['config_file'],
                             '-i',str(Vars['iteration']+1),
                             '-a',str(meanRA),
                             '-d',str(meanDEC),
                             '-z',str(RA_band),
                             '-c', str(DEC_band),
                             '-f', str(meanFreq),
                             '-b',str(F_band),
                             '-T', str(Vars['param_time']),
                             '-s', str(meanF1dot),
                             '-m', str(F_dot_band),
                             '--F2', str(0),
                             '--F2Band', str(0),
                             '--F3', str(0),
                             '--F3Band', str(0),
                             '-t',str(Vars['new_coherence_time'])])
  
  os.system(command_string)

sys.exit(0)
