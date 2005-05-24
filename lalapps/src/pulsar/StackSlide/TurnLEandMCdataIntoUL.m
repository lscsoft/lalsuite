function TurnLEandMCdataIntoUL(jobID,fileList,injectedH0,confidence,numSFT,printToStdOut,outputFile,graphOption)
% Usage: TurnLEandMCdataIntoUL(jobID,fileList,injectedH0,confidence,numSFT,printToStdOut,outputFile,graphOption)
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
% confidence: confidence wanted in ouput upper limit (UL).
% numSFT: number of SFTs used in the search.
% printToStdOut: print info to stdOut
% outputFile: name of output file
% graphOption: if > 0 then display plots
% Finds UL with specified confidence based on LE; histograms the power from the MC Simulation.

% convert parameters when needed:
if (ischar(injectedH0))
    injectedH0=str2num(injectedH0);
end
if (ischar(confidence))
    confidence=str2num(confidence);
end
if (ischar(numSFT))
    numSFT=str2num(numSFT);
end
if (ischar(printToStdOut))
    printToStdOut=str2num(printToStdOut);
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

% Find the search and Monte Carlo (mc) files:
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

%searchFileList = textRead(listSearchXMLfiles,'%s');
searchFileListLength = length(searchFileList);
start_freqs = [];
bands = [];
freqs = [];
RAs = [];
DECs = [];
fDerivs = [];
loudestEventPowers = [];
widths = [];
fileName = char(searchFileList(1)); % get the IFO name from the first file
tbl = readMeta(fileName,'searchsummary_stackslidesfts');
IFO =char(tbl.ifo);
for i=1:searchFileListLength;
  fileName = char(searchFileList(i));
  tbl = readMeta(fileName,'searchsummary_stackslidesfts');  
  start_freqs = [ start_freqs; tbl.start_freq ];
  bands = [ bands; tbl.band ];
  tbl = readMeta(fileName,'sngl_stackslideperiodic');
  freqs = [ freqs; tbl.frequency ];
  RAs = [ RAs; tbl.sky_ra ];
  DECs = [ DECs; tbl.sky_dec ];
  fDerivs = [ fDerivs; tbl.fderiv_1 ];
  loudestEventPowers = [loudestEventPowers; tbl.power];
  widths = [ widths; tbl.width_bins ];    
end

% Basic info about the search
startFreq = min(start_freqs);
[maxFreq, iMaxFreq] = max(start_freqs);
searchBand = maxFreq + bands(iMaxFreq) - startFreq;
 
% Find overall loudest event:
[loudestEvent, iLE] = max(loudestEventPowers);

% Read in the Monte Carlo Simulation Files:
%mcFileList = textRead(listMCxmlFiles,'%s');
mcFileListLength = length(mcFileList);
power = [];
for i=1:mcFileListLength;
  fileName = char(mcFileList(i));
  tbl = readMeta(fileName,'sngl_stackslideperiodic');
  power = [ power; tbl.power ];
end
powerLength = length(power);

% find Loudest Event
power = sort(power);
% Which bin corresponds to falseDismissalRate = 1 - confidence:
iFDRateBin = floor((1.0-confidence)*powerLength + 0.5); % nearest bin
iFDRateBin = iFDRateBin + 1; % go one bin above this and adjust this to be LE.
% Find h0 that would make this equal to loudest event; this corresponds
% to UL at specified confidence, with uncertainty of 1/sqrt(numSFTs).
d2 = 2.0*(power(iFDRateBin) - 1.0);
h0UpperLimit = sqrt( 2.0*(loudestEvent - 1.0)*injectedH0*injectedH0/d2 );
% find uncertainty in Loudest Event
percentUncertainty = (1.0/sqrt(numSFT))/d2;
uncertainty = percentUncertainty*h0UpperLimit; 

% Open file for output
fid = fopen(outputFile,'w');

if (printToStdOut > 0)
 fprintf('\nIFO, startFreq, band = %s, %20.10f, %20.10f\n\n',IFO, startFreq, searchBand);
end
fprintf(fid,'\nIFO, startFreq, band = %s, %20.10f, %20.10f\n\n',IFO, startFreq, searchBand);

if (printToStdOut > 0)
 fprintf('Loudest event: RA, DEC, fDeriv1, f, power, width\n');
 fprintf('%20.10f %20.10f %24.10e %20.10f %20.10f %12.1f\n\n',RAs(iLE),DECs(iLE),fDerivs(iLE),freqs(iLE),loudestEventPowers(iLE),widths(iLE));
end 
fprintf(fid,'Loudest event: RA, DEC, fDeriv1, f, power, width\n');
fprintf(fid,'%20.10f %20.10f %24.10e %20.10f %20.10f %12.1f\n\n',RAs(iLE),DECs(iLE),fDerivs(iLE),freqs(iLE),loudestEventPowers(iLE),widths(iLE));

if (printToStdOut > 0)
 fprintf('The estimated %g confidence Upper Limit = %g +/- %g \n\n',100.0*confidence,h0UpperLimit,uncertainty);
end
fprintf(fid,'The estimated %g confidence Upper Limit = %g +/- %g \n\n',100.0*confidence,h0UpperLimit,uncertainty);

if (printToStdOut > 0)
 fprintf('Loudest Event, Start Frequency, Search Band, Estimated Upper Limit, Uncertainty = \n');
 fprintf('%20.10f %20.10f %20.10f %20.10e %20.10e\n',loudestEventPowers(iLE),startFreq,searchBand,h0UpperLimit,uncertainty);
end
fprintf(fid,'Loudest Event, Start Frequency, Search Band, Estimated Upper Limit, Uncertainty = \n');
fprintf(fid,'%20.10f %20.10f %20.10f %20.10e %20.10e\n',loudestEventPowers(iLE),startFreq,searchBand,h0UpperLimit,uncertainty);

if graphOption > 0
  % Adjust power to h0UpperLimit
  power = (power - 1.0)*h0UpperLimit*h0UpperLimit/(injectedH0*injectedH0);
  power = power + 1.0;
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
