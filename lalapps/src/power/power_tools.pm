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
package power_tools;

#make %STATUS available to calling scripts
use Exporter; 
@ISA = qw(Exporter); 
@EXPORT = qw(	%STATUS 
					f_parseParametersFile 
					f_setupOutputDirs 
					f_updateNotebook
					f_updateNotebookPage
					f_getRunNumber
					f_getDateYYMMDD
					f_buildCacheFile); 
				
my %STATUS = (
					P => "Pending",
					R => "Running",
					C => "Complete",
					U => "User Review Required",
					E => "Error",
					BC => "Bad Cache File",
					NF => "Output file not found");
										
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
#   f_updateNotebook 
#-----------------------------------------------------------------------------------
#  - "keeps the notebook " by creating a html page that lists the
#		runnumbers and links to pages with more details on specific runs
#-----------------------------------------------------------------------------------
#  Returns null
#-----------------------------------------------------------------------------------
sub f_updateNotebook {
	my ($params,$path) = @_;
	
	open RUN_TABLE, "${$params}{'OUTPUT_PATH'}/" . ${$params}{'RUN_TABLE'} 
			or die "Couln't open path/$params{'JOBS_TABLE'}\n";
	
	my %runTable = ();
	
	#load the table into a hash array for sorting, use runid as key and user,desc as value
	#  ??? is used as an arbitrary delimitter
	while(<RUN_TABLE>){
		my ($id,$uid, $desc) = split "\t";
		$runTable{$id} = "$uid\t$desc";
	}
	close RUN_TABLE;
	
	my $tableRows;
	foreach (reverse sort keys %runTable){
		my ($uid, $desc) = split "\t", $runTable{$_};
		$tableRows .= "\t<tr><td nobreak=\"nobreak\"><a href=\"$_/index.html\">$_</a></td><td>$uid</td><td>$desc</td></tr>\n";
	}
	my $notebook = ${$params}{'HTML_PATH'} . ${$params}{'NOTEBOOK'};
	open NOTEBOOK, ">$notebook" or die "Couldn't open $notebook.\n";
	
	print NOTEBOOK << "HTML_PAGE";
<h1>Power Code Notebook</h1>	
<table border="1" cellpadding="2">
<tr><td nobreak="nobreak"><b>Run ID</b></td><td><b>User ID</b></td><td><b>Description</b></td></tr>
$tableRows
</table>
HTML_PAGE

	close NOTEBOOK;
}

#-----------------------------------------------------------------------------------
#   f_updateNotebookPage 
#-----------------------------------------------------------------------------------
#  - Create an HTML page with information on parameters and number of records
#     that is specific to this run
#-----------------------------------------------------------------------------------
#  Returns null
#-----------------------------------------------------------------------------------
sub f_updateNotebookPage {

	if(scalar(@_) != 3){die "Wrong # arguments passed to f_udateNotebookPage.\n";}
	my ($params,$path, $runNum) = @_;
	
	open RUN_TABLE, "$path/${$params}{'JOBS_TABLE'}" 
			or die "Couln't open path/${$params}{'JOBS_TABLE'}\n";
	my $i=0;
	my %statusCount = ( );
	my ($amountOfData, $amountProcessed) = (0,0);
	while(<RUN_TABLE>){
		chomp;
		my ($statusCode, $statusDescription, $startSec, $endSec, $cache, $xml) = split "\t";
		$i++;
		if ($statusCount{$statusDescription}){
			$statusCount{$statusDescription}++;
		}else{
			$statusCount{$statusDescription} = 1;
		}
		$amountOfData += $endSec - $startSec;
		if($statusCode eq "C"){$amountProcessed += $endSec - $startSec;}
	}
	
	my $statusRows;
	foreach(sort keys %statusCount){
		$statusRows .= "\t<tr><td nobreak=\"nobreak\">$_</td><td>$statusCount{$_}</td></tr>\n";
	}

	my $parametersRows;
	foreach (sort keys %{$params}){
		$tableRows .= "\t<tr><td><b>$_</b></td><td>${$params}{$_}</td></tr>\n";
	}
	
	#create the folder for the html page
	my $htmlPath = ${$params}{'HTML_PATH'} . "/$runNum"; 
	if(! -d $htmlPath) { mkdir $htmlPath or die "Couldn't create $htmlPath.\n";}
	
	my $notebookPage = "$htmlPath/index.html";
	open NOTEBOOK_PAGE, ">$notebookPage" or die "Couldn't open $notebookPage.\n";
	
	print NOTEBOOK_PAGE << "HTML_PAGE";
<h1>$runNum Summary</h1>
<table border="1" cellpadding="2">
	<tr><td><b>Total Number of Jobs</b></td><td>$i</td></tr>
	$statusRows
	<tr><td><b>Number of Seconds in Run</b></td><td>$amountOfData</td></tr>
	<tr><td><b>Number of Seconds Completed</b></td><td>$amountProcessed</td></tr>
</table>
<br /><br />
<table border="1" cellpadding="2">
<tr><td colspan="2" align="center"><b>PARAMETERS</b></td></tr>
$tableRows
</table>
HTML_PAGE

	close NOTEBOOK;
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

#-----------------------------------------------------------------------------------
#   f_buildCacheFile
#-----------------------------------------------------------------------------------
#  - Checks to see if the cache file $framechache exists and has data. If it 
#	doesn't exist or doesn't have data, then function makes system call to 
#   LALdataFind to build the cachefile. Sometimes LALdataFind creates empty cachefiles
#   an additional check is made to make sure the new cachefile has data.
# 
#-----------------------------------------------------------------------------------
#  Returns the size in bytes of the cachefile
#-----------------------------------------------------------------------------------
sub f_buildCacheFile {
	my ($startSec, $stopSec, $framecache, $instrument) = @_;
	
	#print "start=$startSec\nstop=$stopSec\ninstrument=$instrument\nFrameCache=$framecache\n";
	#only call LALdataFind if the cache file doesn't currently exist.
	unless(-f $framecache and -s $framecache > 0){
	
		#error END OF FRAME data sometimes occurs for unknown reasons. Making the
		# cache file 1 second longer fixes the problem
		my $cmd =  "LALdataFind --lal-cache --instrument $instrument --type RDS_R_L1 " .
		 		" --start $startSec --end " . ($stopSec + 1) . " > $framecache";
		print "$cmd\n";
		system $cmd;	
		
		unless (-f $framecache && -s $framecache != 0) {
			print "Warning: Framecache file $framecache was not created. Either the data " .
				"was not found or LALdataFind failed.\n";
			unlink ($framecache);
			return 0;
		}
	}
	return (-s $framecache);
}

#return a true value to idicate the end of the module
return 1;