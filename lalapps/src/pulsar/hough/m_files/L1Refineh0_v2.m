% $Id$
%only valid for initial 10 h_0 values
% 
% script that look into L1kkoutput and create a file with the refine H0 values

clear

Detector = 'H1';
fileinput = strcat(Detector,'kkoutput');
fileoutput = strcat(Detector,'h0newband');
fid = fopen(fileoutput, 'w');


%original h0 values all the same


h0file = strcat('/local_data/badkri/S2-clean/MC_allsky/', Detector);
h0file = strcat(h0file, '_coarse/MC_');
h0file = strcat(h0file, Detector);
h0file = strcat(h0file, '_0_h0');
h0val=load( h0file );

CofH = load(fileinput);
kk =size(CofH);

   nbands = kk(1);
   nh = kk(2)-4;

fmin = CofH(:,2);
fmax = CofH(:,3);
Nmax = CofH(:,4);

for j=1:nbands
  Ch0= CofH(j, 5: kk(2));
  h0min = h0val(nh-1);
  h0max = 5.0*h0val(nh);
  
  if( Ch0(nh)> 0.955)
   small = find( Ch0<0.945);
   h0min = h0val(length(small));
   large = find(Ch0 > 0.955);
   h0max = h0val(large(1));
  end
%  fprintf(fid,'%d %d\n', h0min, h0max );
  fprintf(fid,'%d %d %d %d %d %d\n', j-1, fmin(j), fmax(j), Nmax(j), h0min, h0max );

end

fclose(fid);
