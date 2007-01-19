function plotSpecAvgOutput(inputFileName,outputFileName,chanName,tStart,tEnd,fStart,fEnd,effTBase,deltaFTicks,taveFlag)
% usage: plotSpecAvgOutput(inputFileName,outputFileName,chanName,tStart,tEnd,fStart,fEnd,effTBase,deltaFTicks,taveFlag)
% 
% inputFileName  -- the name of the file with data to input; this file is the output from spec_avg.
% outputFileName -- the name of the output plot to generate (this code adds .pdf and .png extensions and outputs these)
% chanName       -- the name of the channel used to generate the spectrogram data
% tStart         -- the start time of the data
% tEnd           -- the end time of the data
% fStart         -- the start frequency of the data
% fEnd           -- the end frequency of the data
% effTBase       -- the effective time base line: 1/effTBase gives the frequency resolution of the plot
%                   so that a effTBase of 10 seconds means deltaF = 0.1 Hz.
% deltaFTicks    -- the change in frequency between major tick marks (e.g., 5 Hz)
% taveFlag       -- if > 0 then produce StackSlide style time average output without sliding.
%                   The spectrograms is given as subplot(2,1,1) and this spectrum as subplot(2,1,2).
%                   Also, .txt is appended to the plot name, and f versus the normalized averaged
%                   power is output into this text file.

% Convert relevant strings to numbers.
if (ischar(tStart))
    tStart=str2num(tStart);
end
if (ischar(tEnd))
    tEnd=str2num(tEnd);
end
if (ischar(fStart))
    fStart=str2num(fStart);
end
if (ischar(fEnd))
    fEnd=str2num(fEnd);
end
if (ischar(effTBase))
    effTBase=str2num(effTBase);
end
if (ischar(deltaFTicks))
    deltaFTicks=str2num(deltaFTicks);
end
if (ischar(taveFlag))
    taveFlag=str2num(taveFlag);
end

xIn = load(inputFileName);

x = transpose(log10(xIn));

figure(1);
if (taveFlag > 0)
   subplot(2,1,1)
end
imagesc(flipud(x));
% Show ticks every deltaFTicks Hz:
deltaTicks = deltaFTicks*effTBase;
numTicks = floor((fEnd - fStart)/deltaFTicks);
vecTicks = 1 + deltaTicks*(0:numTicks);
vecTicks(numTicks + 1) = vecTicks(numTicks + 1) - 1; % for purpose of plot, adjust last tick 
vecFLabels = fStart + deltaFTicks*(0:numTicks);
%set(gca, 'YTick', [ 1 51 101 151 201 251 301 351 401 451 500 ]);
%set(gca, 'YTickLabel', fliplr([50 55 60 65 70 75 80 85 90 95 100 ]));
set(gca, 'YTick', vecTicks);
set(gca, 'YTickLabel', fliplr(vecFLabels));
titleString = sprintf('Spectrogram for %s; GPS %d - %d s.',chanName,tStart,tEnd);
title(titleString);
xlabel('SFT number (see table for corresponding date)');
ylabel('Frequency (Hz)');
colorbar;

if (taveFlag > 0)
  subplot(2,1,2)
  % Produce StackSlide style time average output without sliding:
  timeaverageFileName = sprintf('%s_timeaverage',inputFileName);
  [fk, xout] = textread(timeaverageFileName,'%f %f');
  semilogy(fk,xout)
  %set(gca, 'XTick', vecTicks);
  %set(gca, 'XTickLabel', vecFLabels);
  titleString = sprintf('Spectrum for %s; averaged over GPS %d - %d s.',chanName,tStart,tEnd);
  title(titleString);
  ylabel('Normalized Average Power');
  xlabel('Frequency (Hz)');
  outputTextFile = sprintf('%s.txt',outputFileName);
  fid = fopen(outputTextFile,'w');
  kMax = length(xout);
  for k = 1:kMax
      fprintf(fid,'%f %f\n',fk(k),xout(k));
  end
  fclose(fid);
end
outputPDFFileName = sprintf('%s.pdf',outputFileName);
saveas(1, outputPDFFileName);
outputPNGFileName = sprintf('%s.png',outputFileName);
saveas(1, outputPNGFileName);
delete(1);

return;
