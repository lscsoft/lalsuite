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
% effTBaseFull   -- 1/effTBaseFull gives the frequency resolution of the average normalized spectra plots. This value is the timbaseline of the sfts.
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


xIn = load(filename);%cg; loads the contents of the file into Xin, this should be a matlab binary file according to matlab help, but hte file is ascii!?!?!?!?  Does not seam to matter.  This line loads my test1 data for 3 SFTs inot Xin no probs.  Typing in on the command line you need to put the filename in 'filename'.  Xin is then a 3 by 100 array of doubles.
y = flipud(transpose((xIn)));%cg; flipud flipes the array/matrix upside down, i.e. reverses the row direction of the array.  Transpose obvisouly transposes the array.  So the array gets put up on its side and then turned upside down.  This is so the values appear in the correct way when displaying the array as an image, which is what imagesc does.

%cg;these lines get the relvant info out of the file name
%---------------
undrscr = findstr('_',filename);                                 
fStart = str2num(filename((undrscr(1)+1):(undrscr(2)-1)));      % start frequency 
fEnd = str2num(filename((undrscr(2)+1):(undrscr(3)-1)));        % end frequency
tStart = str2num(filename((undrscr(4)+1):(undrscr(5)-1)));      % start time
tEnd = str2num(filename((undrscr(5)+1):end));                   % end time
%----------------

% calculate characteristic values   cg;  This bit of code gets the data from the input file and orders it appropriately into a matrix for the imagesc matlab function.

y_temp = [ ];
for ii=1:length(y(:,1)); %cg; starts the for loop, starting value for ii is 1, end value for ii is length(y(:,1)).  The length function returns the length of the vector in question, or if it is an array the length of its longest dimension.
  y_temp = [y_temp,y(ii,:)];
end



cutoffval = median(y_temp)+5*(median(y_temp)/sqrt(effTBase));%cg; this cut off value is used in the loop below to get rid of any outliers.  Not sure about the formula used for this calc.  But be aware, will need to change the value being sqrt'ed.Old value is 180, changed to effTBase
%maximum = max(max(y));
%minimum = min(min(y));
     
% replace every value more than three stddeviations away from mean-value

for ii=1:length(y(:,1));
       for jj=1:length(y(1,:));
           if y(ii,jj)>= cutoffval;%cg; 
               y(ii,jj)= cutoffval;
           end
       end
end

filename_date=sprintf('%s_date',filename);
fid_date = fopen(filename_date);
datein=textscan(fid_date, '%d %d %d %d %d %d %d');%once you have used textscan you are at the bottom of the file.
fclose(fid_date);
%so datein now contains a number of vectors each containing a number of cells.  To access the vetcor datein{1} to datein{7}, and then
%I can used datein{1}(1) to datein{n}(n)


%cg; plot the spectrogram and the normalised average power.
%--------------------------------------------------------
figure(1);  %cg; creates a figure graphics object, this figure contains both plots
'thisisfig1'
%if (taveFlag > 0)
   %subplot(2,1,1) %cg; subplot creates mulitple plots for adding to the figure. (2,1,1) specifies that the figure will be 2 plots deep by one plot across, and this is the first plot.  If the line was subplot(2,1,1);plot(graph), then this would plot graph.  This is only used if taveFlag > 0 and we intend to plot the normalised average power.  Will need to get rid of the subplot bit for my modified code.
%end

imagesc((y));%cg;  this line creates the specgram.  Lines after just label it up.  Plots the matrix y.
xlabel('SFT number (see table for corresponding date)','FontSize',13);
ylabel('frequency [Hz]','FontSize',13);
% Show ticks every deltaFTicks Hz:
deltaFTicks = deltaFTicks/10.00;
deltaTicks = deltaFTicks * (effTBaseFull/effTBase);%cg; leave this line in, note its deltaTicks not deltaFTicks

numTicks = floor((fEnd - fStart)/deltaFTicks);%cg; sorts out the labelling of the frequencies on the y axis
vecTicks = 1 + deltaTicks*(0:numTicks); %cg; creates vecticks,for ftsart = 56, fend= 57. from 0 to numticks (0-10), in 1.0, 1.2, 1.4 ...3.0
vecTicks(numTicks + 1) = vecTicks(numTicks + 1) - 1; % for purpose of plot, adjust last tick. cg; for 1hz range this gets rid of everything from 2.1 to 3.0.
vecFLabels = fStart + deltaFTicks*(0:numTicks); %cg;this takes array 1.0, 1.2 ...2.0, multiplies by deltaFTicks (0.1) for this case, and adds to sfstart, so we now have 56.0, 56.1
set(gca, 'YTick', vecTicks);%cg; gca is current axis handle. think this sets the x axis. vecTicks is 1,2,3,4,5, etc. or is it?
%set(gca, 'YTick', vecTicks);
set(gca, 'YTickLabel', fliplr(vecFLabels));
titleString = sprintf('Spectrogram for %s; GPS %d - %d s.' ,chanName,tStart,tEnd);
title(titleString,'Interpreter','none');
colorbar('eastoutside');

%-------------------------------------------
%cg, now print the outputs to a file.
print('-dpng',[outputFileName '.png']);
set(gcf, 'PaperOrientation', 'landscape');%cg; gcf is current figure handle.
set(gcf, 'PaperPosition', [0 0 11 8.5]);
print('-dpdf','-loose',[outputFileName '.pdf'])

delete(1);%closes the figures.

if (taveFlag > 0) %cg; the code below only executes if taveFlag > 0, this bit plots the normalised average power
  %subplot(2,1,2)cg; if printing the specgram and graph seperately, dont need this line.
  figure(2);
%cg; this bit creates the other two text files, these both have the suffix .txt on the end.
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
      if xout(jj)>=cutoff;%cg; xout is one of the properties plotted.
          xout(jj)=cutoff;
      end
  end
 

%cg; Plot the normalised average power vs freq
%-------------------------------------------------
%cg; this uses the blah_blah_blah_timeaverage file output fropm spec_avg.c
  plot(fk,xout)
  titleString = sprintf('Spectrum for %s; averaged over GPS %d - %d s.',chanName,tStart,tEnd);
  title(titleString,'Interpreter','none');
  ylabel('Normalized Average Power');
  xlabel('Frequency (Hz)');
  
  fclose(fid);
  fclose(fid2);

  %cd; these lines put in to plot the normalised avergae power spectrum into a seperate file to the specgram.
  %---------------



  %--------------

  if (thresholdSNR > 0)
  % input the reference file and look for coincidence lines above thresholdSNR
  %----------------------------------------------------------------------------
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
         if (jMax > lengthSNRout); jMax = lengthSNRout; end;
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
  %---------------------------------------------------------------------------------
  
  outputFileName2 = sprintf('%s_2',outputFileName);
	
  print('-dpng',[outputFileName2 '.png'])
  set(gcf, 'PaperOrientation', 'landscape');
  set(gcf, 'PaperPosition', [0 0 11 8.5]);
  print('-dpdf','-loose',[outputFileName2 '.pdf'])
  delete(2);
end

%-------------------------------------------
%cg; this bit just does the final outputs, the png file and the pdf file of the spectrogram

%  print('-dpng',[outputFileName '.png'])
%  set(gcf, 'PaperOrientation', 'landscape');
%  set(gcf, 'PaperPosition', [0 0 11 8.5]);
%  print('-dpdf','-loose',[outputFileName '.pdf'])






return;
