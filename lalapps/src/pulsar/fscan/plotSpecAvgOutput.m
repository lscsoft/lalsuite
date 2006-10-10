function plotSpecAvgOutput(inputFileName,outputFileName,chanName,tStart,tEnd,fStart,fEnd,effTBase,deltaFTicks)
% usage: plotSpecAvgOutput(inputFileName,outputFileName,chanName,tStart,tEnd,fStart,fEnd,effTBase,deltaFTicks)
%

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

xIn = load(inputFileName);

x = transpose(log10(xIn));

figure(1);
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
saveas(1, outputFileName);
