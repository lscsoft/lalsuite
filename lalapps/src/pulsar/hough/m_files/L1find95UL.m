% $Id$
%
%  Matlab script to get  approx h0(95%) by linear interpolation
%
% Remember to add the path where this file is located:
% addpath /scratch/sintes/CVSDIR/waves/people/sintes/PULSAR/CODES/m_files/
%   To run simply type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

Detector = 'H1';

fileoutput = strcat(Detector,'95UL');
fid = fopen(fileoutput, 'w');

fileinput = strcat(Detector,'outCH2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getting the different bands, loudest event per band
% and  h0 and C(h0)

CofH = load(fileinput);
kk =size(CofH);

   nbands = kk(1);
   nh = (kk(2)-4)/2;

fmin = CofH(:,2);
fmax = CofH(:,3);
Nmax = CofH(:,4);

for j=1:nbands
  Ch0= CofH(j, 5: kk(2));
  CHnew = reshape(Ch0, 2,nh);
  h0vec = CHnew(1,:);
  CLvec = CHnew(2,:);
 
  UL=0;
  h01=h0vec(1);
  CL1=CLvec(1);
 
  indices = find(CLvec > 0.95);
  if( length(indices)>0)
    k=indices(1);
    h02=h0vec(k);
    CL2=CLvec(k);
    if(k>1)
      h01=h0vec(k-1);
      CL1=CLvec(k-1);
    else
      h01 = h0vec(1);
      CL1 = CLvec(1);
      h02 = h0vec(2);
      CL2 = CLvec(2);
    end
      slope = (h02-h01)/(CL2-CL1);
      UL=h01 +slope*(0.95 -CL1);
  else
    indices = find(CLvec < 0.95);
    k = length(indices);
    h02 = h0vec(k);
    CL2 = CLvec(k);
    h01 = h0vec(k-1);
    CL1 = CLvec(k-1);  
    slope = (h02 - h01)/(CL2 -CL1);
    UL = h01 + slope * (0.95 - CL1);
  end

   

 fprintf(fid,'%d %d %d %d %d \n', j-1, fmin(j), fmax(j), Nmax(j), UL );
end

fclose(fid);
