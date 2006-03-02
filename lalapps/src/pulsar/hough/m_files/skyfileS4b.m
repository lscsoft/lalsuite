% $Id$
%
%  Matlab script to create one skyfile 

fileoutput= 'skyfileS4c';
fid = fopen(fileoutput, 'w');

%north pole
 fprintf(fid,'%7.5f %7.5f %7.5f %7.5f \n', 0.0000, pi/2.0, 0.4, 0.4 );
 
%delta 1.179 \Delta alpha = 1.0472 6 patches
 
delta = 1.179;
DRa=1.0472;
 
for j=0:5
  fprintf(fid,'%7.5f %7.5f %7.5f %7.5f \n', j*DRa, delta, DRa*cos(delta), 0.393 );
end

%delta 0.786 \Delta alpha = 0.5236 , 12 patches
 
delta = 0.786;
DRa=0.5236;
 
for j=0:11
  fprintf(fid,'%7.5f %7.5f %7.5f %7.5f \n', j*DRa, delta, DRa*cos(delta), 0.393 );
end

%delta 0.393 \Delta alpha = 0.393 , 16 patches
 
delta = 0.393;
DRa=0.393;
 
for j=0:15
  fprintf(fid,'%7.5f %7.5f %7.5f %7.5f \n', j*DRa, delta, DRa*cos(delta), 0.393 );
end


%delta 0.0 \Delta alpha = 0.393 , 16 patches
 
delta = 0.0;
DRa=0.393;
 
for j=0:15
  fprintf(fid,'%7.5f %7.5f %7.5f %7.5f \n', j*DRa, delta, DRa*cos(delta), 0.393 );
end


%delta -0.393 \Delta alpha = 0.393 , 16 patches
 
delta = -0.393;
DRa=0.393;
 
for j=0:15
  fprintf(fid,'%7.5f %7.5f %7.5f %7.5f \n', j*DRa, delta, DRa*cos(delta), 0.393 );
end


%delta -0.786 \Delta alpha = 0.5236 , 12 patches
 
delta = -0.786;
DRa=0.5236;
 
for j=0:11
  fprintf(fid,'%7.5f %7.5f %7.5f %7.5f \n', j*DRa, delta, DRa*cos(delta), 0.393 );
end

%delta -1.179 \Delta alpha = 1.0472 6 patches
 
delta = -1.179;
DRa=1.0472;
 
for j=0:5
  fprintf(fid,'%7.5f %7.5f %7.5f %7.5f \n', j*DRa, delta, DRa*cos(delta), 0.393 );
end


%south pole
 fprintf(fid,'%7.5f %7.5f %7.5f %7.5f \n', 0.0000, -pi/2.0, 0.4, 0.4 );
 
fclose(fid);

