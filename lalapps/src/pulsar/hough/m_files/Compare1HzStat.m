
% $Id$
%
%  Matlab script to plot statistics and compare them with theoretical 
%  distributions
%
% Remember to add the path where this file is located:
% addpath /scratch/sintes/CVSDIR/waves/people/sintes/PULSAR/CODES/m_files/
%   To run simply type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

Directory = '/scratch/sintes/S2/Driver_allsky/';
%Detector = 'L1';
Detector = 'H1';
DirName = '/skypatch_';
prefix ='/freq_';
sufix = '_histo';
sufix2 = '_stats';


skypatch = 9;
bandnumber =2;  % the band of interest


OutDir = strcat(Directory,Detector);
ThisDir = strcat(DirName, int2str( skypatch ));
ThisDir = strcat(OutDir,ThisDir);
file = strcat(prefix, int2str( bandnumber ));

     file1 = strcat(file, sufix);
     file1 = strcat(ThisDir,file1);

     file2 = strcat(file, sufix2);
     file2 = strcat(ThisDir,file2);

alpha = exp(-1.6)

     output=load(file1);
     Nsft = length(output)-1;
     partialSum = sum( output(:,2));
     x = output(:,1);
     y = binopdf(x,Nsft, alpha);
   bar(x, output(:,2)/partialSum)
   hold on
   plot(x,y, 'ro')
   
   %legend('binopdf(n,671,0.20)','p(n|h=0) experiment')
   %title('L1 patch 9, 202-203 Hz' )
   legend('binopdf(n,1886,0.20)','p(n|h=0) experiment')
   title('H1 patch 9, 202-203 Hz' )
  
   NC=load(file2);
   Nmax= NC(:,4);
   Nmin= NC(:,5);
   Nav = NC(:,6);
   Nstd= NC(:,7);
   Nfre= NC(:,8);
   
   figure
    
   plot(Nfre,Nmax,'b+')
   hold on
    plot(Nfre,Nav,'g+')
    plot(Nfre,Nmin,'r+')
    plot(Nfre,Nstd,'m+')
 % title('L1 patch 9, 202-203 Hz' )
 title('H1 patch 9, 202-203 Hz' )
