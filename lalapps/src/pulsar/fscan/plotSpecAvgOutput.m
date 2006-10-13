function plotSpecAvgOutput(inputFileName,outputFileName,chanName,tStart,tEnd,fStart,fEnd,effTBase,deltaFTicks,medBins)
% usage: plotSpecAvgOutput(inputFileName,outputFileName,chanName,tStart,tEnd,fStart,fEnd,effTBase,deltaFTicks,medBins)
% 
% inputFileName  -- the name of the file with data to input; this file is the output from spec_avg.
% outputFileName -- the name of the output plot to generate (the extension determines the type, e.g., .pdf)
% chanName       -- the name of the channel used to generate the spectrogram data
% tStart         -- the start time of the data
% tEnd           -- the end time of the data
% fStart         -- the start frequency of the data
% fEnd           -- the end frequency of the data
% effTBase       -- the effective time base line: 1/effTBase gives the frequency resolution of the plot
%                   so that a effTBase of 10 seconds means deltaF = 0.1 Hz.
% deltaFTicks    -- the change in frequency between major tick marks (e.g., 5 Hz)
% medBins        -- if > 0 then this many frequency bins is used to find a running median normalized averaged
%                   StackSlide style output spectrum. The spectrograms is given as subplot(2,1,1) and this
%                   spectrum as subplot(2,1,2).  Also, .txt is appended to the plot name, and f versus the
%                   normalized averaged power is output into this text file.

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
if (ischar(medBins))
    medBins=str2num(medBins);
end

xIn = load(inputFileName);

x = transpose(log10(xIn));

figure(1);
if (medBins > 0)
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

if (medBins > 0)
  subplot(2,1,2)
  % Produce StackSlide style output:
  % normalize each row using a median value every 10 bins:
  jMax = length(xIn(:,1));
  kMax = length(xIn(1,:));
  medBinso2 = floor(medBins/2);
  kMaxmmedBinso2 = kMax - medBinso2;
  kMaxmmedBins = kMax - medBins;
  for j = 1:jMax
    for k = 1:kMax
        if ( (k > medBinso2) & (k < kMaxmmedBinso2) )
           xss(j,k) = xIn(j,k)/ median( xIn(j,k - medBinso2:k + medBinso2) );
        elseif ( (k <= medBinso2) & (k < kMaxmmedBins) )
           xss(j,k) = xIn(j,k)/ median( xIn(j,k:k + medBins) );
        elseif ( (k > medBins) & (k >= kMaxmmedBinso2) )
           xss(j,k) = xIn(j,k)/ median( xIn(j,k - medBins:k) );
        else
           xss(j,k) = xIn(j,k)/ median( xIn(j,:) );
        end
    end
  end
  for k = 1:kMax
    xout(k) = sum(xss(:,k))/kMax;
  end
  fk = fStart + (0:kMax-1)/effTBase;
  semilogy(fk,xout)
  %set(gca, 'XTick', vecTicks);
  %set(gca, 'XTickLabel', vecFLabels);
  titleString = sprintf('Spectrum for %s; averaged over GPS %d - %d s.',chanName,tStart,tEnd);
  title(titleString);
  ylabel('Normalized Average Power');
  xlabel('Frequency (Hz)');
  outputTextFile = sprintf('%s.txt',outputFileName);
  fid = fopen(outputTextFile,'w');  
  for k = 1:kMax
      fprintf(fid,'%f %f\n',fk(k),xout(k));
  end
  fclose(fid);
end
saveas(1, outputFileName);

return;
