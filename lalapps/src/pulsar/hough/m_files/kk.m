
fileoutput = 'Hough_MultiIFO_UL.txt';
fid = fopen(fileoutput, 'w');

Nbands = length(fmin);

for bandnumber = 1:Nbands
  fprintf(fid,'%d  %d \n', fmin(bandnumber), UL(bandnumber) );
end
fclose(fid);
