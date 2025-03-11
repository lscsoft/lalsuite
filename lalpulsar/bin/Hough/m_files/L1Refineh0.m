%only valind for initial 10 h_0 values
%
% script that look into L1kkoutput and create a file with the refine H0 values

clear

Detector = 'L1';
fileinput = strcat(Detector,'kkoutput');
fileoutput = strcat(Detector,'h0newband');
fid = fopen(fileoutput, 'w');


%original h0 values

h0val=(1:10)*1e-23;

CofH = load(fileinput);
kk =size(CofH);

   nbands = kk(1);
   nh = kk(2)-4;

fmin = CofH(:,2);
fmax = CofH(:,3);
Nmax = CofH(:,4);

for j=1:nbands
  Ch0= CofH(j, 5: kk(2));
  h0min = h0val(9);
  h0max = 5.0e-22;

  if( Ch0(nh)> 0.955)
   small = find( Ch0<0.945);
   h0min = h0val(length(small));
   large = find(Ch0 > 0.955);
   h0max = h0val(large(1));
  end
   fprintf(fid,'%d %d %d %d %d %d\n', j-1, fmin(j), fmax(j), ...
           Nmax(j), h0min, h0max );
end

fclose(fid);
