%
% $Id$
%
%  Matlab script to plot C(h_0) for the 1st iteration
%
% Remember to add the path where this file is located:
% addpath /scratch/sintes/CVSDIR/waves/people/sintes/PULSAR/CODES/m_files/
%   To run simply type  L1PlotCofH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

Detector = 'L1';
fileinput = strcat(Detector,'kkoutput');

CofH = load(fileinput);
kk =size(CofH);
nh = kk(2)-4;

fmin = CofH(:,2);
fmax = CofH(:,3);
Nmax = CofH(:,4);
nn= length(fmin);

f=[fmin,fmax];
f2=reshape(f', 1,nn*2);

figure
 colores=['g', 'y','k','c', 'm', 'b', 'r', 'g', 'y','k','c', 'm', 'b', 'r', 'g', 'b','k'];
 
 
for j=5:kk(2)
  cc = CofH(:,j);
  ccc= [cc,cc];
  c2 =  reshape(ccc', 1,nn*2);
   plot(f2, c2,colores(j) )
  %plot(f2, c2,'-','LineWidth',j/9)
  hold on
end
axis([200 400 -0.05 1.05]);
xlabel ('Frequency [Hz]'); 
ylabel ('C(h_0)'); 
title('S2: L1 5.0e-23 to 1.4e-22')
%legend('All sky', 'North belt')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

