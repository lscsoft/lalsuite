% $Id$
%
%  Matlab script to plot approx h0(95%) by linear interpolation
%
% Remember to add the path where this file is located:
% addpath /scratch/sintes/CVSDIR/waves/people/sintes/PULSAR/CODES/m_files/
%   To run simply type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

Detector = 'L1';
fileinput = strcat(Detector,'95UL');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


UL =  load(fileinput);

fmin = UL(:,2);
fmax = UL(:,3);
Nmax = UL(:,4);
 
upper = UL(:,5);
 
Nbands =length(fmin);

totalband = 0;

for bn=1:Nbands
  plot( [fmin(bn) fmax(bn)], [upper(bn) upper(bn)]);
  hold on
totalband(bn) = fmax(bn) - fmin(bn);
end

  sum(totalband) - 0.02*6*200

plot((fmin+fmax)/2, upper, '+')
 grid
 xlabel ('Frequency [Hz]');
 ylabel ('h_0 95%');
 titlestring = strcat('S2-', Detector);
 titlestring = strcat(titlestring, ': Preliminary Upper Limit');
 title(titlestring)
