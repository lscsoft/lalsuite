#-----------------------------------------------------------------------------------
# power_tools.pm
#-----------------------------------------------------------------------------------
#  This module contains the perl functions shared by
#  by the power scripts which initially consisted of 
#  createJobsTable, processJobsTable, and countJobsInTable
# 
# Revision History
# July, 2003 - Dennis Mackin <dsmackin@stupidlinux.com>  - Created module
# $Id$
#-----------------------------------------------------------------------------------
					
					
#-----------------------------------------------------------------------------------
#  f_parseParametersFile
#-----------------------------------------------------------------------------------
#  - parses the parameters file
#  - dies if parameters file is invalid
#-----------------------------------------------------------------------------------
#  RETURNS : the address of the parameters hash array 
#-----------------------------------------------------------------------------------
sub f_parseParametersFile{

	my $parametersFile  = shift ;
	
	open PARAMS, $parametersFile or die "In  f_parseParametersFile: Couldn't open $parametersFile.";
	my %params = ();
	
	while(<PARAMS>){
		#print $_,"\n";
		#get rid of whitespace from beginning and end of string
		s/^\s*//;
		s/$\s*//;
		#print $_,"\n";		
		# skip comments
		if ( /^#/ ){ next;} 
		
		#remove the newline character
		chomp;

		#if length is 0 move to next line
		if (length($_) < 1){next;}
		
		my @pair  = split "=";
		
		#if a parameter/value pair is not found, then exit the program
		if(scalar(@pair) != 2) {die "failed to find parameter value pair in f_parseParametersFile .\n";}
		
		#clean off any " marks
		$pair[1] =~ s/"//g;
		
		#clean off whitespace characters from the parameter name
		$pair[0] =~ s/\s//g;
		$pair[1] =~ s/\s*$//g;
		$pair[1] =~ s/^\s*//g;
		
		$params{$pair[0]} = $pair[1];
	}
	
	close PARAMS;
	return \%params;
}

#-----------------------------------------------------------------------------------
#   f_setupOutputDirs()
#-----------------------------------------------------------------------------------
#  - checks to see if output dirs exist. if not, creates them
#-----------------------------------------------------------------------------------
#  Returns 
#-----------------------------------------------------------------------------------
sub f_setupOutputDirs {
	my ($path,$date,$runNum) = @_;

	$path =   "$path/$date-$runNum"; 
	my $xmldir = "$path/xml/";
	my $logdir = "$path/log/";
	my $errdir = "$path/err/";
	my $outdir = "$path/out/";
	
	if(! -d $path) { mkdir $path or die "Couldn't create $path.\n";}
	if(! -d $xmldir) { mkdir $xmldir or die "Couldn't create $xmldir.\n";}
	if(! -d $logdir) { mkdir $logdir or die "Couldn't create $logdir\n";}
	if(! -d $errdir) { mkdir $errdir or die "Couldn't create $errdir\n";}
	if(! -d $outdir) { mkdir $outdir or die "Couldn't create $outdir\n";}			
	print "built paths . . .\n";
}


#-----------------------------------------------------------------------------------
#  f_getRunNumber
#-----------------------------------------------------------------------------------
#  - Looks in the run output directory to determine the next run
#    number for the current day
#-----------------------------------------------------------------------------------
#  Returns the next run number as NNN
#-----------------------------------------------------------------------------------
sub f_getRunNumber {
	my ($path, $date) = @_;
	opendir PATH, $path	;
	
	#get a bottom to top list of the files that end with DATE-XXX
	my $searchString = $date . '-\d{3}';
	my @dirs = reverse sort grep /$\$searchString/, readdir PATH;
	closedir PATH;
	
	if(scalar @dirs > 0){
		my @parts  = reverse split "-", $dirs[0];
		my $maxNum = $parts[0];
		my $nextNum = $maxNum + 1;
		
		#make sure nextNum is like NNN, not NN or N
		while(length($nextNum) < 3){ $nextNum = "0" . $nextNum;}
		
		return $nextNum;
	}

	#by default return 001
	return "001";
}


#-----------------------------------------------------------------------------------
#  f_getDateYYMMDD
#-----------------------------------------------------------------------------------
#  Returns date as YYMMDD
#-----------------------------------------------------------------------------------
sub f_getDateYYMMDD{
    my %month = (
                Jan => '01',
                Feb => '02',
                Mar => '03',
                Apr => '04',
                May => '05',
                Jun => '06',
                Jul => '07',
                Aug => '08',
                Sep => '09',
                Oct => '10',
                Nov => '11',
                Dec => '12'
                );
    my @timeparts = split(" ", localtime());

    #add in leading 0 for day if needed
    if ($timeparts[2]<10) {$timeparts[2] = "0" . $timeparts[2];};

	 #change YYYY to YY
	 $timeparts[4] =~ s/^20//;

    return "$timeparts[4]$month{$timeparts[1]}$timeparts[2]";
}

#return a true value to idicate the end of the module
return 1;