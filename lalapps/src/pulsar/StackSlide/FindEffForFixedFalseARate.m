function FindEffForFixedFalseARate(listSearchXMLfiles,listMCxmlFiles,injectedH0,falseAlarmRate,numSFT,inputThreshold,outputFile,graphOption)
% Usage: FindEffForFixedFalseARate(listSearchXMLfiles,listMCxmlFiles,injectedH0,falseAlarmRate,numSFT,inputThreshold,outputFile,graphOption)
% searchXMLfile: file with list of files with output of StackSlide Search
% listMCxmlFiles: file with list of files with output of corresponding Monte Carlo Simulation
% injectedH0: value of signal amplitude h_0 injected during Monte Carlo Simulation
% falseAlarmRate: falseAlarmRate used to determine threshold for efficiency measurement..
% numSFT: Not used.
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
if (ischar(numSFT))
    numSFT=str2num(numSFT);
end
if (ischar(inputThreshold))
    inputThreshold=str2num(inputThreshold);
end
if (ischar(graphOption))
    graphOption=str2num(graphOption);
end
% Read in the search results files:
searchFileList = textRead(listSearchXMLfiles,'%s');
searchFileListLength = length(searchFileList);
start_freqs = [];
bands = [];
freqs = [];
RAs = [];
DECs = [];
fDerivs = [];
eventPowers = [];
widths = [];
fileName = char(searchFileList(1)); % get the IFO name for the first file
tbl = readMeta(fileName,'searchsummary_stackslidesfts');
IFO =char(tbl.ifo);
totalNumTemplates = 0.0;
for i=1:searchFileListLength;
  fileName = char(searchFileList(i));
  tbl = readMeta(fileName,'searchsummary_stackslidesfts');  
  start_freqs = [ start_freqs; tbl.start_freq ];
  bands = [ bands; tbl.band ];
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

% Open file for output
fid = fopen(outputFile,'w');

% Print basic info about the search
startFreq = min(start_freqs);
[maxFreq, iMaxFreq] = max(start_freqs);
searchBand = maxFreq + bands(iMaxFreq) - startFreq;
fprintf('\nIFO, startFreq, band = %s, %g, %g\n',IFO, startFreq,searchBand);
fprintf(fid,'\nIFO, startFreq, band = %s, %g, %g\n',IFO, startFreq,searchBand);
 
% Set up overall loudest event, and print loudest event for each frequency band.
if (inputThreshold > 0)
  falseARateThreshold = inputThreshold;
  fprintf('\nfalseARateThreshold was input as = %g\n\n',falseARateThreshold);
  fprintf(fid,'\nfalseARateThreshold was input as = %g\n\n',falseARateThreshold);
else 
 eventPowers = sort(eventPowers);
 lengthEventPowers = length(eventPowers);
 iFARateBin = lengthEventPowers - (falseAlarmRate*totalNumTemplates)
 falseARateThreshold = eventPowers(iFARateBin)
 fprintf('\nfalseARateThreshold found is: %g \n',falseARateThreshold);
 fprintf(fid,'\nfalseARateThreshold event found is: %g \n',falseARateThreshold);
end

% Read in the Monte Carlo Simulation Files:
fileList = textRead(listMCxmlFiles,'%s');
fileListLength = length(fileList);
power = [];
for i=1:fileListLength;
  fileName = char(fileList(i));
  tbl = readMeta(fileName,'sngl_stackslideperiodic');
  power = [ power; tbl.power ];
end
powerLength = length(power);

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
  %plot(x,y,'+');
  %C = (numSTK^numSTK)/factorial(numSTK-1);
  %C = 1;
  %xplot = min(x) + [0:960]*(max(x) - min(x))/999.0;
  %yPDF = C*xplot.^(numSFT-1).*exp(-1.0*numSFT*xplot);
  %yPDFint = trapz(xplot,yPDF);
  %yint = trapz(x,y);
  %yPDF = yint*yPDF/yPDFint;
  %xlabel('STACKSLIDE POWER');
  %ylabel('NUMBER');
  %hold
  %plot(xplot,yPDF,'r');
  %hold
end

fclose(fid);
return;
