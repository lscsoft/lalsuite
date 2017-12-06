function ReportOnULsFromAtoBCondorDag(commentString,fileList,outputFile,confidence,outputOption,printToStdOut,graphOption,xAccOverc,yAccOverc,zAccOverc)
% Usage: ReportOnULsFromAtoBCondorDag(commentString,fileList,outputFile,confidence,outputOption,printToStdOut,graphOption,xAccOverc,yAccOverc,zAccOverc)
% commentString: comment pertaining to this search (e.g, S4 run)
% fileLIst: entire list of files for all jobs.  This files has the format <node letter+jobID> <filename>, e.g,:
%
%           A0 filenameA0_0.txt
%
%           B0 filenameB0_0.txt
%
%           A1 filenameA1_0.txt
%
%           B1 filenameB1_0.txt
%
%           etc...
%
% mcFileList is found by matching files in fileList that are B jobs.
% confidence: desired confidence of ULs.
% outputFile: name of output file (give [] if no output file)
% if (bitand(outputOption,1) > 0) will compute the S and SBins values for the loudest event.
% if (bitand(outputOption,2) > 0)
%     When reading in Monte Carlo results to find the linear least squares fit to h_0 vs C,
%     start at j=2 below; this skip the first row of the searchresults_stackslidemontecarlo
%     which is redundant under some options. Note that if the first row of the
%     searchresults_stackslidemontecarlo is estimated UL it will be skipped in any case.
% if (bitand(outputOption,4) > 0) and (graphOption > 0) show estimated UL and conf on ULs vs Conf graph.
% graphOption: if > 0 then display plots, including a plot of ULs vs Confidence for band(i) for i=graphOption

if (ischar(confidence))
    confidence=str2num(confidence);
end
if (ischar(outputOption))
    outputOption=str2num(outputOption);
end
if (ischar(printToStdOut))
    printToStdOut=str2num(printToStdOut);
end
if (ischar(graphOption))
    graphOption=str2num(graphOption);
end

% set up jobIDs for A and B jobs:
strA = sprintf('A');
strB = sprintf('B');

% Read in entire list of files:
[jobIDList, entireFileList] = textread(fileList,'%s %s');
entireFileListLength = length(entireFileList);

% Find search files with loudest event (le) from A jobs:
countA = 0;
for i=1:entireFileListLength
    if ( length( strfind(char(jobIDList(i)),strA) ) > 0 )
       countA = countA + 1;
       leFileList{countA} = char(entireFileList(i));
    end
end

% Find Monte Carlo (mc) files from B jobs:
countB = 0;
for i=1:entireFileListLength
    if ( length( strfind(char(jobIDList(i)),strB) ) > 0 )
       countB = countB + 1;
       mcFileList{countB} = char(entireFileList(i));
    end
end

% Read in the Monte Carlo Simulation Files:
mcFileListLength = length(mcFileList);
leFileListLength = length(leFileList);
if (mcFileListLength ~= leFileListLength)
   disp('Error: Number of LE files not equal to number of MC files.')
   return;
end
% get the IFO name from the first file
fileName = char(mcFileList(1));
tbl = readMeta(fileName,'searchsummary_stackslidesfts');
IFO =char(tbl.ifo);
start_time = tbl.start_time;
start_time_ns = tbl.start_time_ns;
duration = tbl.duration;
sft_baseline = tbl.sft_baseline;
num_sfts = tbl.num_sfts;

% initial array with loudest event (loudest power) information
lePWR = [];
leRA = [];
leDEC = [];
leFREQ = [];
leFDOT = [];
leS = [];
leSHz = [];
leSBins = [];

start_freq = [];
band = [];
% First row of searchresults_stackslidemontecarlo contains estimated ULs
upper_limit_est = [];
confidence_est = [];
converged_est = [];
% Search all rows for each Job to find best ULs
best_upper_limit = [];
delta_upper_limit = [];
best_confidence = [];
delta_confidence= [];
best_converged = []; % not really used in this code
for i=1:leFileListLength;
  
  % Get the loudest event
  fileName = char(leFileList(i));
  tbl = readMeta(fileName,'sngl_stackslideperiodic');
  [maxPWR, imax] = max(tbl.power); % Just in case more than one event is in the table , get loudest
  
  lePWR = [lePWR; maxPWR];
  leRA = [leRA; tbl.sky_ra(imax)];
  leDEC = [leDEC; tbl.sky_dec(imax)];
  leFREQ = [leFREQ; tbl.frequency(imax)];
  leFDOT = [leFDOT; tbl.fderiv_1(imax)];

  if (bitand(outputOption,1) > 0)
      [S, SHz, Sbins] = FindSkyPosAndFDotSValue(xAccOverc,yAccOverc,zAccOverc,tbl.sky_ra(imax),tbl.sky_dec(imax),tbl.frequency(imax),tbl.fderiv_1(imax),duration,sft_baseline);
  else
      S = 0.0;
      SHz = 0.0;
      Sbins = 0;
  end    
  leS = [leS; S];
  leSHz = [leSHz;SHz];
  leSBins = [leSBins;round(Sbins)];

  % Get the MC data
  fileName = char(mcFileList(i));
  tbl = readMeta(fileName,'searchresults_stackslidemontecarlo');
  tbl.loudest_event(1);                                     % same for all rows per job
  start_freq = [start_freq; tbl.start_freq(1)];             % same for all rows per job
  band = [band; tbl.band(1)];                               % same for all rows per job
  upper_limit_est = [upper_limit_est; tbl.upper_limit(1)];  % estimates are in the first row
  confidence_est = [confidence_est; tbl.confidence(1)];     % estimates are in the first row
  converged_est = [converged_est; tbl.converged(1)];        % estimates are in the first row    

  [mindiff, imin] = min(  abs(tbl.confidence - confidence)  );   % find row closest to desired confidence
  uls = [];
  confs = [];
  if (bitand(outputOption,2) > 0)
    jstart = 2;
    % When reading in Monte Carlo results to find the linear least squares fit to h_0 vs C,
    % start at j=2 below; this skip the first row of the searchresults_stackslidemontecarlo
    % which is redundant under some options. Note that if the first row of the
    % searchresults_stackslidemontecarlo is estimated UL it will be skipped in any case.
  else
    jstart = 1;
  end
  for j=jstart:length(tbl.upper_limit)
      if ((tbl.confidence(j) > -1.0) && (tbl.converged(j) > -1))
         uls = [uls; tbl.upper_limit(j)];
         confs = [confs; tbl.confidence(j)];
      end
  end
  [P,S] = polyfit(confs,uls,1);
  [thisUL,deltaUL] = polyval(P,confidence,S);
  [P,S] = polyfit(uls/mean(uls),confs,1);  % run the fit with uls and conf reversed to get an error estimate on the conf.
  [Y,deltaConf] = polyval(P,thisUL/mean(uls),S);
  if i==graphOption
    figure(4);
    y = polyval(P,uls/mean(uls),S);
    [ySorted,isort] = sort(y);
    ulsSorted = uls(isort);    
    if (bitand(outputOption,4) > 0)
      plot(uls,confs,'+',ulsSorted,ySorted,'k',thisUL,confidence,'b*',tbl.upper_limit(2),tbl.confidence(2),'o');
      legend('measurements','best fit line','best fit UL','estimated UL');
    else
      plot(uls,confs,'+',ulsSorted,ySorted,'k',thisUL,confidence,'b*');
      legend('measurements','best fit line','best fit UL');
    end
    xlabel('h_0');
    ylabel('confidence');
    titleString = sprintf('StackSlide %s h_0 vs. conf. from MC simulation, and best fit. %s Band: %g-%g Hz.',IFO,commentString,start_freq(i),start_freq(i)+band(i));
    title(titleString);
  end  
  best_upper_limit = [best_upper_limit; thisUL];  
  delta_upper_limit = [delta_upper_limit; deltaUL];
  best_confidence = [best_confidence; confidence];
  delta_confidence = [delta_confidence; deltaConf];  
  best_converged = [best_converged; -1];
end

% Open file for output
if length(outputFile > 0)
  fid = fopen(outputFile,'w');
  outputFileNbrs = strcat(outputFile,'.nbrs');
  fidNbrs = fopen(outputFileNbrs,'w');
else
  fid =  0
end

if (printToStdOut > 0)
 fprintf('StackSlide %s %4.1f percent conf. upper limits for %9.4f to %9.4f Hz band. %s\n',IFO,100.0*confidence,min(start_freq),max(start_freq),commentString);

 fprintf('\nStart freq  Loudest event  Estimated UL\n\n');
 for i=1:mcFileListLength
   fprintf('%15.4f %15.4f %20.10e\n',start_freq(i),lePWR(i),upper_limit_est(i));
 end
  
 fprintf('\n                    Loudest Event Parameters                                              Upper Limit        Confidence');
 fprintf('\n Start_freq    RA      DEC         FREQ        FDOT             S       SBins  PWR        UL    +/-    Unc    Conf  +/-  Unc\n\n');
 for i=1:mcFileListLength
   fprintf('%10.4f %9.6f %9.6f %12.6f %14.6e %12.4e %4i %7.4f   %8.2e +/- %8.2e  %5.2f +/- %6.3f\n',start_freq(i),leRA(i),leDEC(i),leFREQ(i),leFDOT(i),leS(i),leSBins(i),lePWR(i),best_upper_limit(i),delta_upper_limit(i),best_confidence(i),delta_confidence(i));
 end
 
 [maxdiff, imax] = max( abs(best_upper_limit - upper_limit_est)./best_upper_limit );
 fprintf('\nMax fractional difference between ULs and estimated ULs is %15.6f for %15.4f Hz\n',maxdiff,start_freq(imax));

 [minUL, imin] = min( best_upper_limit );
 fprintf('Best UL is %20.10e for %15.4f Hz\n',minUL,start_freq(imin));

 [maxUL, imax] = max( best_upper_limit );
 fprintf('Worst ULs is %20.10e for %15.4f Hz\n',maxUL,start_freq(imax)); 
end

if (fid > 0)
 fprintf(fid,'StackSlide %s %4.1f percent conf. upper limits for %9.4f to %9.4f Hz band. %s\n',IFO,100.0*confidence,min(start_freq),max(start_freq),commentString);
 
 fprintf(fid,'\nStart freq  Loudest event  Estimated UL\n\n');
 for i=1:mcFileListLength
   fprintf(fid,'%15.4f %15.4f %20.10e\n',start_freq(i),lePWR(i),upper_limit_est(i));
 end
  
 fprintf(fid,'\n                    Loudest Event Parameters                                              Upper Limit        Confidence');
 fprintf(fid,'\n Start_freq    RA      DEC         FREQ        FDOT             S       SBins  PWR        UL    +/-    Unc    Conf  +/-  Unc\n\n');
 for i=1:mcFileListLength
   fprintf(fid,'%10.4f %9.6f %9.6f %12.6f %14.6e %12.4e %4i %7.4f   %8.2e +/- %8.2e  %5.2f +/- %6.3f\n',start_freq(i),leRA(i),leDEC(i),leFREQ(i),leFDOT(i),leS(i),leSBins(i),lePWR(i),best_upper_limit(i),delta_upper_limit(i),best_confidence(i),delta_confidence(i));
 end 
 % Just print out the number for later use
 for i=1:mcFileListLength
   fprintf(fidNbrs,'%10.4f %9.6f %9.6f %12.6f %14.6e %12.4e %4i %7.4f   %8.2e  %8.2e  %5.2f  %6.3f\n',start_freq(i),leRA(i),leDEC(i),leFREQ(i),leFDOT(i),leS(i),leSBins(i),lePWR(i),best_upper_limit(i),delta_upper_limit(i),best_confidence(i),delta_confidence(i));
 end

 [maxdiff, imax] = max( abs(best_upper_limit - upper_limit_est)./best_upper_limit );
 fprintf(fid,'\nMax fractional difference between ULs and estimated ULs is %15.6f for %15.4f Hz\n',maxdiff,start_freq(imax));

 [minUL, imin] = min( best_upper_limit );
 fprintf(fid,'Best UL is %20.10e for %15.4f Hz\n',minUL,start_freq(imin));

 [maxUL, imax] = max( best_upper_limit );
 fprintf(fid,'Worst ULs is %20.10e for %15.4f Hz\n',maxUL,start_freq(imax));

 fclose(fid);
 fclose(fidNbrs);
end

if graphOption > 0
   figure(1);
     subplot(2,1,1);
       plot(start_freq,lePWR,'*')
       xlabel('Frequency (Hz)');
       ylabel('Loudest StackSlide Power');
       titleString = sprintf('StackSlide %s loudest events. %s',IFO,commentString);
       title(titleString);
     subplot(2,1,2);     
       plot(start_freq,upper_limit_est,'*');
       xlabel('Frequency (Hz)');
       ylabel('Estimated Upper Limit');
       titleString = sprintf('StackSlide %s estimated %4.1f percent conf. ULs. %s',IFO,100.0*confidence,commentString);
       title(titleString);     
   figure(2);
     subplot(2,1,1);
       plot(start_freq,lePWR,'*')
       xlabel('Frequency (Hz)');
       ylabel('Loudest StackSlide Power');
       titleString = sprintf('StackSlide %s loudest events. %s',IFO,commentString);
       title(titleString);
     subplot(2,1,2);
       plot(start_freq,best_upper_limit,'*');
       xlabel('Frequency (Hz)');
       ylabel('Upper Limit');
       titleString = sprintf('StackSlide %s %4.1f percent conf. upper limits. %s',IFO,100.0*confidence,commentString);
       title(titleString);       
    figure(3);
       PlotRaDecPwrAllSkyLEs(leRA,leDEC,lePWR,1,50,1);
end

return;
