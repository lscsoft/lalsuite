function ReportOnULsFromIteratedMC(run,skyBand,fileList,outputFile,confidence,outputOption,printToStdOut,graphOption)
% Usage: ReportOnULsFromIteratedMC(fileList,outputFile,confidence,outputOption,printToStdOut,graphOption)
% run: which run (e.g., 'S3')
% skyBand: which band on the sky (e.g., '0')
% fileLIst: entire list of files for all jobs.  This files has the format <node letter+jobID> <filename>, e.g,:
%
%           A0 filenameA0_0.txt
%           B0 filenameB0_0.txt
%           B0 filenameB0_1.txt
%           C0 filenameC0_0.txt
%           D0 filenameD0_0.txt
%
%           A1 filenameA1_0.txt
%           B1 filenameB1_0.txt
%           B1 filenameB1_1.txt
%           C1 filenameC1_0.txt
%           D1 filenameD1_0.txt
%
%           etc...
%
% mcFileList is found by matching files in fileList that are D jobs.
% confidence: desired confidence of ULs
% outputFile: name of output file (give [] if no output file)
% outputOption: what to output
% graphOption: if > 0 then display plots

if (ischar(confidence))
    outputOption=str2num(confidence);
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
strD = sprintf('D');

% Read in entire list of files:
[jobIDList, entireFileList] = textread(fileList,'%s %s');
entireFileListLength = length(entireFileList);

% Find Monte Carlo (mc) files from D jobs:
countD = 0;
for i=1:entireFileListLength
    if ( length( strfind(char(jobIDList(i)),strD) ) > 0 )
       countD = countD + 1;
       mcFileList{countD} = char(entireFileList(i));
    end
end

fileName = char(mcFileList(1)); % get the IFO name from the first file
tbl = readMeta(fileName,'searchsummary_stackslidesfts');
IFO =char(tbl.ifo);

% Read in the Monte Carlo Simulation Files:
mcFileListLength = length(mcFileList);
loudest_event = [];
start_freq = [];
%band = [];
% First row of searchresults_stackslidemontecarlo contains estimated ULs
upper_limit_est = [];
confidence_est = [];
converged_est = [];
% Search all rows for each Job to find best ULs
best_upper_limit = [];
best_confidence = [];
best_converged = [];
for i=1:mcFileListLength;
  fileName = char(mcFileList(i));
  tbl = readMeta(fileName,'searchresults_stackslidemontecarlo');
  loudest_event = [loudest_event; tbl.loudest_event(1)];    % same for all rows per job
  start_freq = [start_freq; tbl.start_freq(1)];             % same for all rows per job
  upper_limit_est = [upper_limit_est; tbl.upper_limit(1)];  % estimates are in the first row
  confidence_est = [confidence_est; tbl.confidence(1)];     % estimates are in the first row
  converged_est = [converged_est; tbl.converged(1)];        % estimates are in the first row    

  [mindiff, imin] = min(  abs(tbl.confidence - confidence)  );   % find row closest to desired confidence
  
  best_upper_limit = [best_upper_limit; tbl.upper_limit(imin)];  % from row closest to desired confidence
  best_confidence = [best_confidence; tbl.confidence(imin)]; % from row closest to desired confidence
  best_converged = [best_converged; tbl.converged(imin)];    % from row closest to desired confidence
end

% Open file for output
if length(outputFile > 0)
  fid = fopen(outputFile,'w');
else
  fid =  0
end

if (printToStdOut > 0)
 fprintf('Upper limits for StackSlide %s %s data, for %15.4f to %15.4f Hz band and sky band %s.\n',run,IFO,min(start_freq),max(start_freq),skyBand);

 fprintf('\nStart freq  Loudest event  Estimated UL  Confidence  Converged \n\n');
 for i=1:mcFileListLength
   fprintf('%15.4f %15.4f %20.10e %15.6f %i\n',start_freq(i),loudest_event(i),upper_limit_est(i),confidence_est(i),converged_est(i));
 end
 percent_converged =(1.0*length(find(converged_est>0)))/(1.0*mcFileListLength);
 fprintf('Percent of estimated ULs that converged to desired confidence = %15.6f\n',percent_converged);
 [maxdiff, imax] = max(  abs(confidence_est - confidence)  );
 fprintf('Max divergence from desired confidence for estimated ULs is %15.6f for %15.4f Hz\n',maxdiff,start_freq(imax)); 
 
 fprintf('\nStart freq  Loudest event  UL  Confidence  Converged \n\n');
 for i=1:mcFileListLength
   fprintf('%15.4f %15.4f %20.10e %15.6f %i\n',start_freq(i),loudest_event(i),best_upper_limit(i),best_confidence(i),best_converged(i));
 end
 percent_converged =(1.0*length(find(best_converged>0)))/(1.0*mcFileListLength);
 fprintf('Percent of ULs that converged to desired confidence = %15.6f\n',percent_converged);
 [maxdiff, imax] = max(  abs(best_confidence - confidence)  );
 fprintf('Max divergence from desired confidence for ULs is %15.6f for %15.4f Hz\n',maxdiff,start_freq(imax));
 
 [maxdiff, imax] = max( abs(best_upper_limit - upper_limit_est)./best_upper_limit );
 fprintf('\nMax fractional difference between ULs and estimated ULs is %15.6f for %15.4f Hz\n',maxdiff,start_freq(imax));

 [minUL, imin] = min( best_upper_limit );
 fprintf('Best UL is %20.10e for %15.4f Hz\n',minUL,start_freq(imin));

 [maxUL, imax] = max( best_upper_limit );
 fprintf('Worst ULs is %20.10e for %15.4f Hz\n',maxUL,start_freq(imax));
 
end

if (fid > 0)
 fprintf(fid,'Upper limits for StackSlide %s %s data, for %15.4f to %15.4f Hz band and sky band %s.\n',run,IFO,min(start_freq),max(start_freq),skyBand);
 
 fprintf(fid,'\nStart freq  Loudest event  Estimated UL  Confidence  Converged \n\n');
 for i=1:mcFileListLength
   fprintf(fid,'%15.4f %15.4f %20.10e %15.6f %i\n',start_freq(i),loudest_event(i),upper_limit_est(i),confidence_est(i),converged_est(i));
 end
 percent_converged =(1.0*length(find(converged_est>0)))/(1.0*mcFileListLength);
 fprintf(fid,'Percent of estimated ULs that converged to desired confidence = %15.6f\n',percent_converged);
 [maxdiff, imax] = max(  abs(confidence_est - confidence)  );
 fprintf(fid,'Max divergence from desired confidence for estimated ULs is %15.6f for %15.4f Hz\n',maxdiff,start_freq(imax)); 
 
 fprintf('fid,\nStart freq  Loudest event  UL  Confidence  Converged \n\n');
 for i=1:mcFileListLength
   fprintf(fid,'%15.4f %15.4f %20.10e %15.6f %i\n',start_freq(i),loudest_event(i),best_upper_limit(i),best_confidence(i),best_converged(i));
 end
 percent_converged =(1.0*length(find(best_converged>0)))/(1.0*mcFileListLength);
 fprintf(fid,'Percent of ULs that converged to desired confidence = %15.6f\n',percent_converged);
 [maxdiff, imax] = max(  abs(best_confidence - confidence)  );
 fprintf(fid,'Max divergence from desired confidence for ULs is %15.6f for %15.4f Hz\n',maxdiff,start_freq(imax));
 
 [maxdiff, imax] = max( abs(best_upper_limit - upper_limit_est)./best_upper_limit );
 fprintf(fid,'\nMax fractional difference between ULs and estimated ULs is %15.6f for %15.4f Hz\n',maxdiff,start_freq(imax));

 [minUL, imin] = min( best_upper_limit );
 fprintf(fid,'Best UL is %20.10e for %15.4f Hz\n',minUL,start_freq(imin));

 [maxUL, imax] = max( best_upper_limit );
 fprintf(fid,'Worst ULs is %20.10e for %15.4f Hz\n',maxUL,start_freq(imax));

 fclose(fid);
end

if graphOption > 0
   figure(1);
     plot(start_freq,loudest_event,'*')
     xlabel('Frequency (Hz)');
     ylabel('Loudest StackSlide Power');
     titleString = sprintf('StackSlide %s %s Loudest Events For Sky Band %s',run,IFO,skyBand);
     title(titleString);
   figure(2);   
     subplot(2,1,1);
       plot(start_freq,upper_limit_est,'*');
       xlabel('Frequency (Hz)');
       ylabel('Estimated Upper Limit');
       titleString = sprintf('StackSlide %s %s Estimated ULs For Sky Band %s',run,IFO,skyBand);
       title(titleString);
     subplot(2,1,2);   
       plot(start_freq,confidence_est,'*');
       xlabel('Frequency (Hz)');
       ylabel('Confidence of Estimated UL');
       titleString = sprintf('StackSlide %s %s Confidence of Estimated ULs For Sky Band %s',run,IFO,skyBand);
       title(titleString);       
   figure(3);
     subplot(2,1,1);
       plot(start_freq,best_upper_limit,'*');
       xlabel('Frequency (Hz)');
       ylabel('Upper Limit');
       titleString = sprintf('StackSlide %s %s ULs For Sky Band %s',run,IFO,skyBand);
       title(titleString);       
     subplot(2,1,2);
       plot(start_freq,best_confidence,'*');
       xlabel('Frequency (Hz)');
       ylabel('Confidence');
       titleString = sprintf('StackSlide %s %s Confidence of ULs For Sky Band %s',run,IFO,skyBand);
       title(titleString);              
end

return;
