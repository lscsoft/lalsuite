function TurnLEandMCdataIntoUL(listSearchXMLfiles,listMCxmlFiles,injectedH0,confidence,numSFT,inputLE,outputFile,graphOption)
% Usage: TurnLEandMCdataIntoUL(listSearchXMLfiles,listMCxmlFiles,injectedH0,confidence,numSFT,inputLE)
% searchXMLfile: file with list of files with output of StackSlide Search
% listMCxmlFiles: file with list of files with output of corresponding Monte Carlo Simulation
% injectedH0: value of signal amplitude h_0 injected during Monte Carlo Simulation
% confidence: confidence wanted in ouput upper limit (UL).
% numSFT: number of SFTs used in the search.
% inputLE: if > 0 then uses this as the Loudest Event (LE); otherwise get LE from searchXMLfile.
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
if (ischar(inputLE))
    inputLE=str2num(inputLE);
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
loudestEventPowers = [];
widths = [];
fileName = char(searchFileList(1)); % get the IFO name for the first file
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

% Open file for output
fid = fopen(outputFile,'w');

% Print basic info about the search
startFreq = min(start_freqs);
[maxFreq, iMaxFreq] = max(start_freqs);
searchBand = maxFreq + bands(iMaxFreq) - startFreq;
fprintf('\nIFO, startFreq, band = %s, %g, %g\n',IFO, startFreq,searchBand);
fprintf(fid,'\nIFO, startFreq, band = %s, %g, %g\n',IFO, startFreq,searchBand);
 
%fprintf('\n Table Loudest Events: \n\n');
%fprintf(fid,'\n Table Loudest Events: \n\n');
%fprintf('        RA          DEC       fDeriv        frequency       power         width \n');
%fprintf(fid,'        RA          DEC       fDeriv        frequency       power         width \n');
%for i=1:length(freqs)
%  if (loudestEventPowers(i) > 0.0)
%    fprintf('%12.6f %12.6f %15.6e %12.6f %12.6f %12.6f\n',RAs(i),DECs(i),fDerivs(i),freqs(i),loudestEventPowers(i),widths(i));
%    fprintf(fid,'%12.6f %12.6f %15.6e %12.6f %12.6f %12.6f\n',RAs(i),DECs(i),fDerivs(i),freqs(i),loudestEventPowers(i),widths(i));
%  end
%end
% Set up overall loudest event, and print loudest event for each frequency band.
if (inputLE > 0)
  loudestEvent = inputLE;
  fprintf('\nLoudest event was input as = %g\n\n',loudestEvent);
  fprintf(fid,'\nLoudest event was input as = %g\n\n',loudestEvent);
else 
 [loudestEvent, i] = max(loudestEventPowers);
 fprintf('\nLoudest event found is: %g for RA, DEC, fDerivs, f, power, width: \n',loudestEvent);
 fprintf(fid,'\nLoudest event found is: %g for RA, DEC, fDerivs, f, power: \n',loudestEvent);
 fprintf('%12.6f %12.6f %15.6e %12.6f %12.6f %12.6f\n\n',RAs(i),DECs(i),fDerivs(i),freqs(i),loudestEventPowers(i),widths(i));
 fprintf(fid,'%12.6f %12.6f %15.6e %12.6f %12.6f %12.6f\n\n',RAs(i),DECs(i),fDerivs(i),freqs(i),loudestEventPowers(i),widths(i));
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

% find Loudest Event
power = sort(power);
% Which bin corresponds to falseDismissalRate = 1 - confidence:
iFDRateBin = floor((1.0-confidence)*powerLength + 0.5); % nearest bin
iFDRateBin = iFDRateBin + 1; % go one bin above this and adjust this to be LE.
% Find h0 that would make this = to loudest event; this corresponds
% to UL at specified confidence, with uncertainty of 1/sqrt(numSFTs).
d2 = 2.0*(power(iFDRateBin) - 1.0)
h0UpperLimit = sqrt( 2.0*(loudestEvent - 1.0)*injectedH0*injectedH0/d2 );
% find uncertainty in Loudest Event
percentUncertainty = (1.0/sqrt(numSFT))/d2;
uncertainty = percentUncertainty*h0UpperLimit; 
fprintf('\nThe %g confidence Upper Limit = %g +/- %g \n',100.0*confidence,h0UpperLimit,uncertainty);
fprintf(fid,'\nThe %g confidence Upper Limit = %g +/- %g \n',100.0*confidence,h0UpperLimit,uncertainty);

fprintf('Frequency and Upper Limit = \n');
fprintf('%20.10f %20.10e\n',startFreq, h0UpperLimit);
fprintf(fid,'Frequency and Upper Limit = \n');
fprintf(fid,'%20.10f %20.10e\n',startFreq, h0UpperLimit);

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
