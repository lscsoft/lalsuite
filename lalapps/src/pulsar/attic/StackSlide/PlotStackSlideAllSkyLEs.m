function PlotStackSlideAllSkyLEs(fileNameOrFileList,nameOrList,decimation,parea,type)
% Usage: PlotStackSlideAllSkyLEs(fileNameOrFileList,nameOrList,decimation,parea,type)
% fileNameOrFileList: An xml file with loudest event vs sky postion or a txt file with a list of such files.
% nameOrList:  if 0 fileNameOrFileList is an xml file, else fileNameOrFileList is a file with a list of such files.
%              if nameOrList == 3  a map of patch colors is shown instead.
% decimation: factor to reduce amount of data
% (e.g., 1 = no reduction; 10 = reduce by factor of 10).
% parea: point area (in pts squared) of dots used in plot.
% type: type of plot (1 = scatter plot, 2 = use projected ra

ra = [];
dec = [];
deriv1 = [];
freq = [];
pwr = [];
if (nameOrList == 0)
  % read in data
  fileName = char(fileNameOrFileList)
  tbl = readMeta(fileName,'sngl_stackslideperiodic');  
  ra = [ ra; tbl.sky_ra ];
  dec = [ dec; tbl.sky_dec ];
  deriv1 = [ deriv1; tbl.fderiv_1 ];
  freq = [ freq; tbl.frequency ];
  pwr = [ pwr; tbl.power ];
else
  % Read in list of Files:
  fileList = textread(fileNameOrFileList,'%s');
  fileListLength = length(fileList);
  for i=1:fileListLength;
    fileName = char(fileList(i));
    tbl = readMeta(fileName,'sngl_stackslideperiodic');
    ra = [ ra; tbl.sky_ra ];
    dec = [ dec; tbl.sky_dec ];
    deriv1 = [ deriv1; tbl.fderiv_1 ];
    freq = [ freq; tbl.frequency ];
    if (nameOrList == 3)    
      pwr = [ pwr; (i-1)*ones(length(tbl.power),1) ]; % special case!
    elseif (nameOrList == 4)
      pwr = [ pwr; (0.05*(i-1))*tbl.power ]; % special case!
    else
      pwr = [ pwr; tbl.power ];
    end
  end
end

% Find Loudest;
[pwrMax, imax ] = max(pwr);
raMax = ra(imax);
decMax = dec(imax);
deriv1Max = deriv1(imax);
freqMax = freq(imax);
fprintf('The loudest event has RA, DEC, FDERIV_1, FREQ, PWR = \n');
fprintf('%23.10f %23.10f %23.10e %23.10f %23.10f\n',raMax, decMax, deriv1Max, freqMax, pwrMax);

vlength = length(pwr);
ra = ra(1:decimation:vlength);
dec = dec(1:decimation:vlength);
pwr = pwr(1:decimation:vlength);
if (type==1)
   scatter(ra,dec,parea,pwr,'filled');
else
   projra = pi*(1 - cos(dec)) + ra.*cos(dec);
   scatter(projra,dec,parea,pwr,'filled');
end
xlabel('RA (radians)');
ylabel('DEC (radians)');
title('StackSlide Power Loudest Events For Each Sky Position');
return;
