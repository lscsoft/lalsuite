function plotSpecAvgOutput(filename,outputFileName,chanName,effTBase,deltaFTicks,taveFlag,effTBaseFull,thresholdSNR,coinDF,referenceFile)
%
% fileName       -- the name of the file with data to input; this file is the output from spec_avg.
% outputFileName -- base name of the output spectrogram files; will append .pdf and .png to this name.
% chanName       -- the name of the channel used to generate the spectrogram data.
% effTBase       -- the effective time base line: 1/effTBase gives the frequency resolution of the plot.
%                   so that a effTBase of 10 seconds means deltaF = 0.1 Hz.
% deltaFTicks    -- the change in frequency between major tick marks (e.g., 5 Hz).
% taveFlag       -- if > 0 then produce StackSlide style time average output without sliding.
%                   The spectrograms is given as subplot(2,1,1) and this spectrum as subplot(2,1,2).
%                   Also, .txt is appended to the plot name, and f versus the normalized averaged
%                   power is output into this text file.
% effTBaseFull   -- 1/effTBaseFull gives the frequency resolution of the average normalized spectra plots.
% thresholdSNR   -- if > 0 then look for coincident lines with the referenceFile spectra above this threshold.
% coinDF         -- window in frequency to use when looking for coincident lines.
% referenceFile  -- base name of the reference file output by spec_avg; will append _timeaverage to this name.

% Convert relevant strings to numbers.
if (ischar(effTBase))
    effTBase=str2num(effTBase);
end
if (ischar(deltaFTicks))
    deltaFTicks=str2num(deltaFTicks);
end
if (ischar(taveFlag))
    taveFlag=str2num(taveFlag);
end
if (ischar(effTBaseFull))
    effTBaseFull=str2num(effTBaseFull);
end
if (ischar(thresholdSNR))
    thresholdSNR=str2num(thresholdSNR);
end
if (ischar(coinDF))
    coinDF=str2num(coinDF);
end


xIn = load(filename);
y = flipud(transpose((xIn)));
   
undrscr = findstr('_',filename);                                 
fStart = str2num(filename((undrscr(1)+1):(undrscr(2)-1)));      % start frequency
fEnd = str2num(filename((undrscr(2)+1):(undrscr(3)-1)));        % end frequency
tStart = str2num(filename((undrscr(4)+1):(undrscr(5)-1)));      % start time
tEnd = str2num(filename((undrscr(5)+1):end));                   % end time

% calculate characteristic values   

y_temp = [ ];
for ii=1:length(y(:,1));
  y_temp = [y_temp,y(ii,:)];
end
cutoffval = median(y_temp)+5*(median(y_temp)/sqrt(180));
%maximum = max(max(y));
%minimum = min(min(y));
     
% replace every value more than three stddeviations away from mean-value

for ii=1:length(y(:,1));
       for jj=1:length(y(1,:));
           if y(ii,jj)>= cutoffval;
               y(ii,jj)= cutoffval;
           end
       end
end

% plot the spectrogram
figure(1);
if (taveFlag > 0)
   subplot(2,1,1)
end
imagesc((y));
xlabel('SFT number (see table for corresponding date)','FontSize',13)
ylabel('frequency [Hz]','FontSize',13)
% Show ticks every deltaFTicks Hz:
deltaTicks = deltaFTicks*effTBase;
numTicks = floor((fEnd - fStart)/deltaFTicks);
vecTicks = 1 + deltaTicks*(0:numTicks);
vecTicks(numTicks + 1) = vecTicks(numTicks + 1) - 1; % for purpose of plot, adjust last tick 
vecFLabels = fStart + deltaFTicks*(0:numTicks);
set(gca, 'YTick', vecTicks);
set(gca, 'YTickLabel', fliplr(vecFLabels));
titleString = sprintf('Spectrogram for %s; GPS %d - %d s.' ,chanName,tStart,tEnd);
title(titleString,'Interpreter','none');
colorbar('eastoutside');

if (taveFlag > 0)
  subplot(2,1,2)
% Produce StackSlide style time average output without sliding:
  timeaverageFileName = sprintf('%s_timeaverage',filename);
  [fk, xout] = textread(timeaverageFileName,'%f %f');
  outputTextFile = sprintf('%s.txt',filename);
  outputSortedTextFile = sprintf('%s_sorted.txt',outputFileName);
  [xoutSorted,iSorted] = sort(xout,'descend');
  fSorted = fk(iSorted);
  fid = fopen(outputTextFile,'w');
  fid2 = fopen(outputSortedTextFile,'w');
  kMax = length(xout);
  for k = 1:kMax
      fprintf(fid,'%f %f\n',fk(k),xout(k));
      fprintf(fid2,'%f %f\n',fSorted(k),xoutSorted(k));
  end

  stdev_xout = std(xout)
  meanval_xout = mean(xout)
  
  % Read in timestamps file to find the number of SFTs used:
  timestampFileName = sprintf('%s_timestamps',filename);
  [ntmp, ttmp] = textread(timestampFileName,'%f %f');
  % Computed expected 5 sigma cutoff for gaussian noise:
  cutoffmax = 1.0 + 5.0/sqrt(length(ntmp));
  % avoid too small of a maximum cutoff:
  if cutoffmax < 2;
     cutoffmax = 2;
  end
  % compute cutoff from data, but do not exceed cutoffmax:
  cutoff = meanval_xout+(5*stdev_xout);
  if cutoff > cutoffmax;
    cutoff = cutoffmax;
  end
  
  for jj=1:length(xout);
      if xout(jj)>=cutoff;
          xout(jj)=cutoff;
      end
  end
 
  plot(fk,xout)
  titleString = sprintf('Spectrum for %s; averaged over GPS %d - %d s.',chanName,tStart,tEnd);
  title(titleString,'Interpreter','none');
  ylabel('Normalized Average Power');
  xlabel('Frequency (Hz)');
  
  fclose(fid);
  fclose(fid2);
  if (thresholdSNR > 0)
  % input the reference file and look for coincidence lines above thresholdSNR.
   timeaverageFileNameRef = sprintf('%s_timeaverage',referenceFile);
   [fRef, xRef] = textread(timeaverageFileNameRef,'%f %f');
   outputTextFileLines = sprintf('%s_lines.txt',outputFileName);
   fid3 = fopen(outputTextFileLines,'w');
   fprintf(fid3,'\n       COINCIDENT LINES       \n');
   fprintf(fid3,' \n');
   fprintf(fid3,'             INPUT                             REFERENCE          \n');
   fprintf(fid3,' Freq. (Hz)   Power      SNR       Freq. (Hz)   Power      SNR    \n');
   fprintf(fid3,' \n');   
   xoutMean = mean(xout);
   oneOverxoutSTD = 1.0/std(xout);
   xRefMean = mean(xRef);
   oneOverxRefSTD = 1.0/std(xRef);
   SNRout = (xout - xoutMean)*oneOverxoutSTD;
   SNRRef = (xRef - xRefMean)*oneOverxRefSTD;
   lengthSNRout = length(SNRout);
   skip = 0;
   iMax = 0;   
   coincidenceBins = ceil(coinDF*effTBaseFull);
   for j = 1:lengthSNRout
       if (skip == 0)
         jMin = j - coincidenceBins;
         if (jMin < 1); jMin = 1; end;
         jMax = j + coincidenceBins;
        [SNRoutmax, iMaxout] = max(SNRout(jMin:jMax));
         iMaxout = jMin + iMaxout - 1;
         [SNRRefmax, iMaxRef] = max(SNRRef(jMin:jMax));
         iMaxRef = jMin + iMaxRef - 1;
         if ( (SNRoutmax >= thresholdSNR) && (SNRRefmax >= thresholdSNR) )
            skip = 1;
            fprintf(fid3,' %11.6f  %9.4f  %7.2f    %11.6f  %9.4f  %7.2f\n',fk(iMaxout),xout(iMaxout),SNRoutmax,fRef(iMaxRef),xRef(iMaxRef),SNRRefmax);
         end
       else
           if ((j - iMaxRef) > coincidenceBins) && ((j - iMaxout) > coincidenceBins); skip = 0; end;
       end
   end
   fclose(fid3);
  end 
end

print('-dpng',[outputFileName '.png'])
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperPosition', [0 0 11 8.5]);
print('-dpdf','-loose',[outputFileName '.pdf'])


delete(1);


return;
