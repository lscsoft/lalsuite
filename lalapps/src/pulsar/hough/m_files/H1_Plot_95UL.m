% $Id$
%
%  Matlab script to plot approx h0(95%) by linear interpolation
%
% Remember to add the path where this file is located:
% addpath /scratch/sintes/CVSDIR/waves/people/sintes/PULSAR/CODES/m_files/
%   To run simply type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

Detector = 'H1';
fileinput = strcat(Detector,'95UL');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


UL =  load(fileinput);

fmin = UL(:,2);
fmax = UL(:,3);
Nmax = UL(:,4);
 
upper = UL(:,5);
 
Nbands =length(fmin);

for bn=1:Nbands
  plot( [fmin(bn) fmax(bn)], [upper(bn) upper(bn)]);
  hold on
end

plot((fmin+fmax)/2, upper, '+')
 grid
 xlabel ('Frequency [Hz]');
 ylabel ('h_0 95%');
 title('S2-H1: Preliminary Upper Limit')
