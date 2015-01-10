function FindEffForFixedFalseARate(jobID,fileList,injectedH0,falseAlarmRate,inputThreshold,outputFile,graphOption)
% Usage: FindEffForFixedFalseARate(jobID,fileList,injectedH0,falseAlarmRate,inputThreshold,outputFile,graphOption)
% jobID: ID number of the job to be processed.
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
% searchFileList is found by matching files in fileList with A+jobID.
% listMCxmlFiles is found by matching files in fileList with B+jobID.
% searchXMLfile: file with list of files with output of StackSlide Search
% listMCxmlFiles: file with list of files with output of corresponding Monte Carlo Simulation
% injectedH0: value of signal amplitude h_0 injected during Monte Carlo Simulation
% falseAlarmRate: falseAlarmRate used to determine threshold for efficiency measurement.
% inputThreshold: if NOT > 0 get threshold from the falseAlarmRate; otherwise use this as the threshold.
% outputFile: name of output file
% graphOption: if > 0 then display plots
% Finds efficiency for specified falseAlarmRate or inputThreshold; histograms the power from the MC Simulation.

% convert parameters when needed:
if (ischar(injectedH0))
    injectedH0=str2num(injectedH0);
end
if (ischar(falseAlarmRate))
    falseAlarmRate=str2num(falseAlarmRate);
end
if (ischar(inputThreshold))
    inputThreshold=str2num(inputThreshold);
end
if (ischar(graphOption))
    graphOption=str2num(graphOption);
end
if (isnumeric(jobID))
   jobID = num2str(jobID);
end

% set up jobIDs for A and B jobs:
jobIDA = sprintf('A%s',jobID);
jobIDB = sprintf('B%s',jobID);

% Read in entire list of files:
[jobIDList, entireFileList] = textread(fileList,'%s %s');
entireFileListLength = length(entireFileList);

% Find the search and Monte Carlo (mc) files. The search files are not needed if inputThreshold > 0.
countA = 0;
countB = 0;
for i=1:entireFileListLength
    if (strcmp(char(jobIDList(i)),jobIDA))
       countA = countA + 1;
       searchFileList{countA} = char(entireFileList(i));
    elseif (strcmp(char(jobIDList(i)),jobIDB))
       countB = countB + 1;
       mcFileList{countB} = char(entireFileList(i));
    end
end

% The search files are not used if inputThreshold > 0.
if (inputThreshold > 0)
  falseARateThreshold = inputThreshold;
else
  if (countA < 1)
    fprintf('\nError: search files must exist if inputThreshold == 0.\n');
    return;
  end
  searchFileListLength = length(searchFileList);
  freqs = [];
  RAs = [];
  DECs = [];
  fDerivs = [];
  eventPowers = [];
  widths = [];
  totalNumTemplates = 0.0;
  for i=1:searchFileListLength;
    fileName = char(searchFileList(i));
    if (exist(fileName,'file') < 1)
       fprintf('\nError: search files must exist if inputThreshold == 0.\n');
       return;
    end
    tbl = readMeta(fileName,'searchsummary_stackslidesfts');  
    numBins = tbl.num_bins
    numSums = tbl.num_sums
    totalNumTemplates = totalNumTemplates + (numBins*numSums)
    tbl = readMeta(fileName,'sngl_stackslideperiodic');
    freqs = [ freqs; tbl.frequency ];
    RAs = [ RAs; tbl.sky_ra ];
    DECs = [ DECs; tbl.sky_dec ];
    fDerivs = [ fDerivs; tbl.fderiv_1 ];
    eventPowers = [eventPowers; tbl.power];
    widths = [ widths; tbl.width_bins ];    
  end
  eventPowers = sort(eventPowers);
  lengthEventPowers = length(eventPowers);
  iFARateBin = lengthEventPowers - round(falseAlarmRate*totalNumTemplates)
  falseARateThreshold = eventPowers(iFARateBin)
end

% Read in the Monte Carlo Simulation Files:
%mcFileList = textRead(listMCxmlFiles,'%s');
mcFileListLength = length(mcFileList);
power = [];
start_freqs = [];
bands = [];
for i=1:mcFileListLength;
  fileName = char(mcFileList(i));
  tbl = readMeta(fileName,'searchsummary_stackslidesfts');
  IFO =char(tbl.ifo); % ifo should be the same in all files
  start_freqs = [ start_freqs; tbl.start_freq ];
  bands = [ bands; tbl.band ];
  tbl = readMeta(fileName,'sngl_stackslideperiodic');
  power = [ power; tbl.power ];  
end
powerLength = length(power);

% Open file for output
fid = fopen(outputFile,'w');

% Print basic info about the search
startFreq = min(start_freqs);
[maxFreq, iMaxFreq] = max(start_freqs);
searchBand = maxFreq + bands(iMaxFreq) - startFreq;
fprintf('\nIFO, startFreq, band = %s, %g, %g\n',IFO, startFreq,searchBand);
fprintf(fid,'\nIFO, startFreq, band = %s, %g, %g\n',IFO, startFreq,searchBand);

if (inputThreshold > 0)
  fprintf('\nfalseARateThreshold was input as = %g\n\n',falseARateThreshold);
  fprintf(fid,'\nfalseARateThreshold was input as = %g\n\n',falseARateThreshold);
else
  fprintf('\nfalseARateThreshold found is: %g \n',falseARateThreshold);
  fprintf(fid,'\nfalseARateThreshold event found is: %g \n',falseARateThreshold);
end

% find efficiency for finding power with injections above search LE: 
louderPowersLength = length(find(power>falseARateThreshold));
theEfficiency = louderPowersLength/powerLength;
fprintf('\nThe start frequency and efficiency are = \n',theEfficiency);
fprintf('%20.10f %20.10f\n',startFreq, theEfficiency);
fprintf(fid,'\nThe start frequency and efficiency are = \n',theEfficiency);
fprintf(fid,'%20.10f %20.10f\n',startFreq, theEfficiency);

if graphOption > 0
  nBins = 50;
  [y,x] = hist(power,nBins);
  hist(power,nBins);
end

fclose(fid);
return;
