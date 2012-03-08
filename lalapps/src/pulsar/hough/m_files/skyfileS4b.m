%
%  Matlab script to create one skyfile 

fileoutput= 'skyfileS4c';
fid = fopen(fileoutput, 'w');


Ddec=0.393;
%north pole
 fprintf(fid,'%7.5f %7.5f %7.5f %7.5f \n', 0.0000, pi/2.0, 0.4, 0.4 );
 
%delta 1.179 = 3*Ddec \Delta alpha = 0.7854 8 patches
 
delta = 1.179;
DRa=0.7854;
 
for j=0:7
  fprintf(fid,'%7.5f %7.5f %7.5f %7.5f \n', j*DRa, delta, DRa*cos(delta-Ddec/2.0), Ddec );
end

%delta 0.786 \Delta alpha =  0.4834, 13 patches
 
delta = 0.786;
DRa=0.4834;
 
for j=0:12
  fprintf(fid,'%7.5f %7.5f %7.5f %7.5f \n', j*DRa, delta, DRa*cos(delta-Ddec/2.0), Ddec );
end

%delta 0.393 \Delta alpha = 0.393 , 16 patches
 
delta = 0.393;
DRa=0.393;
 
for j=0:15
  fprintf(fid,'%7.5f %7.5f %7.5f %7.5f \n', j*DRa, delta, DRa*cos(delta-Ddec/2.0), Ddec );
end


%delta 0.0 \Delta alpha = 0.393 , 16 patches
 
delta = 0.0;
DRa=0.393;
 
for j=0:15
  fprintf(fid,'%7.5f %7.5f %7.5f %7.5f \n', j*DRa, delta, DRa*cos(delta), Ddec );
end


%delta -0.393 \Delta alpha = 0.393 , 16 patches
 
delta = -0.393;
DRa=0.393;
 
for j=0:15
  fprintf(fid,'%7.5f %7.5f %7.5f %7.5f \n', j*DRa, delta, DRa*cos(delta+Ddec/2.0), Ddec );
end


%delta -0.786 \Delta alpha = 0.4834 , 13 patches
 
delta = -0.786;
DRa=0.4834;
 
for j=0:12
  fprintf(fid,'%7.5f %7.5f %7.5f %7.5f \n', j*DRa, delta, DRa*cos(delta+Ddec/2.0), Ddec );
end

%delta -1.179 \Delta alpha = 0.7854 8 patches
 
delta = -1.179;
DRa=0.7854;
 
for j=0:7
  fprintf(fid,'%7.5f %7.5f %7.5f %7.5f \n', j*DRa, delta, DRa*cos(delta+Ddec/2.0), Ddec );
end


%south pole
 fprintf(fid,'%7.5f %7.5f %7.5f %7.5f \n', 0.0000, -pi/2.0, 0.4, 0.4 );
 
fclose(fid);

