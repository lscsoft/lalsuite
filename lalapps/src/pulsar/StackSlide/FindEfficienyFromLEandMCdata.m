function FindEfficiencyFromLEandMCdata(listSearchXMLfiles,listMCxmlFiles,injectedH0,inputLE,outputFile,graphOption)
% Usage: FindEfficiencyFromLEandMCdata(listSearchXMLfiles,listMCxmlFiles,injectedH0,inputLE)
% searchXMLfile: file with list of files with output of StackSlide Search
% listMCxmlFiles: file with list of files with output of corresponding Monte Carlo Simulation
% injectedH0: value of signal amplitude h_0 injected during Monte Carlo Simulation
% inputLE: if > 0 then uses this as the Loudest Event (LE); otherwise get LE from searchXMLfile.
% outputFileName: Name of file to with result output.
% graphOption: if > 0 then genereate histogram of the power from the MC Simulation.
% Finds efficiency of finding injections > LE

% Read in the search results files:
if (ischar(injectedH0))
    injectedH0=str2num(injectedH0);
end
if (ischar(inputLE))
    inputLE=str2num(inputLE);
end
searchFileList = textread(listSearchXMLfiles,'%s');
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
fprintf('\nIFO, startFreq, band = %s, %g, %g\n',IFO,startFreq,searchBand);
fprintf(fid,'\nIFO, startFreq, band = %s, %g, %g\n',IFO,startFreq,searchBand);
 
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
fileList = textread(listMCxmlFiles,'%s');
fileListLength = length(fileList);
power = [];
for i=1:fileListLength;
  fileName = char(fileList(i));
  tbl = readMeta(fileName,'sngl_stackslideperiodic');
  power = [ power; tbl.power ];
end
powerLength = length(power);

% find efficiency for finding power with injections above search LE: 
louderPowersLength = length(find(power>loudestEvent));
theEfficiency = louderPowersLength/powerLength;
fprintf('\nThe start frequency and efficiency are = \n',theEfficiency);
fprintf('%20.10f %20.10f\n',startFreq, theEfficiency);
fprintf(fid,'\nThe start frequency and efficiency are = \n',theEfficiency);
fprintf(fid,'%20.10f %20.10f\n',startFreq, theEfficiency);

fclose(fid);
return;
