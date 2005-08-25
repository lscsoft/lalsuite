clear

path(path,'/home/siemens/ligotools/matlab')

% Coincident efficiencies

cd /scratch2/xavi/CosmicStrings0

H1ctrigfound0 = readMeta('H1cinjfound.xml','sim_burst');
H2ctrigfound0 = readMeta('H2cinjfound.xml','sim_burst');
L1ctrigfound0 = readMeta('L1cinjfound.xml','sim_burst');

H1cinjmade0 = readMeta('H1cinj.xml','sim_burst');
H2cinjmade0 = readMeta('H2cinj.xml','sim_burst');
L1cinjmade0 = readMeta('L1cinj.xml','sim_burst');

cd /scratch2/xavi/CosmicStrings1

H1ctrigfound1 = readMeta('H1cinjfound.xml','sim_burst');
H2ctrigfound1 = readMeta('H2cinjfound.xml','sim_burst');
L1ctrigfound1 = readMeta('L1cinjfound.xml','sim_burst');

H1cinjmade1 = readMeta('H1cinj.xml','sim_burst');
H2cinjmade1 = readMeta('H2cinj.xml','sim_burst');
L1cinjmade1 = readMeta('L1cinj.xml','sim_burst');

cd /scratch2/xavi/CosmicStrings2

H1ctrigfound2 = readMeta('H1cinjfound.xml','sim_burst');
H2ctrigfound2 = readMeta('H2cinjfound.xml','sim_burst');
L1ctrigfound2 = readMeta('L1cinjfound.xml','sim_burst');

H1cinjmade2 = readMeta('H1cinj.xml','sim_burst');
H2cinjmade2 = readMeta('H2cinj.xml','sim_burst');
L1cinjmade2 = readMeta('L1cinj.xml','sim_burst');

cd /scratch2/xavi/CosmicStrings3

H1ctrigfound3 = readMeta('H1cinjfound.xml','sim_burst');
H2ctrigfound3 = readMeta('H2cinjfound.xml','sim_burst');
L1ctrigfound3 = readMeta('L1cinjfound.xml','sim_burst');

H1cinjmade3 = readMeta('H1cinj.xml','sim_burst');
H2cinjmade3 = readMeta('H2cinj.xml','sim_burst');
L1cinjmade3 = readMeta('L1cinj.xml','sim_burst');

cd /scratch2/xavi/CosmicStrings4

H1ctrigfound4 = readMeta('H1cinjfound.xml','sim_burst');
H2ctrigfound4 = readMeta('H2cinjfound.xml','sim_burst');
L1ctrigfound4 = readMeta('L1cinjfound.xml','sim_burst');

H1cinjmade4 = readMeta('H1cinj.xml','sim_burst');
H2cinjmade4 = readMeta('H2cinj.xml','sim_burst');
L1cinjmade4 = readMeta('L1cinj.xml','sim_burst');

cd /scratch2/xavi/CosmicStrings5

H1ctrigfound5 = readMeta('H1cinjfound.xml','sim_burst');
H2ctrigfound5 = readMeta('H2cinjfound.xml','sim_burst');
L1ctrigfound5 = readMeta('L1cinjfound.xml','sim_burst');

H1cinjmade5 = readMeta('H1cinj.xml','sim_burst');
H2cinjmade5 = readMeta('H2cinj.xml','sim_burst');
L1cinjmade5 = readMeta('L1cinj.xml','sim_burst');

cd /scratch2/xavi/CosmicStrings6

H1ctrigfound6 = readMeta('H1cinjfound.xml','sim_burst');
H2ctrigfound6 = readMeta('H2cinjfound.xml','sim_burst');
L1ctrigfound6 = readMeta('L1cinjfound.xml','sim_burst');

H1cinjmade6 = readMeta('H1cinj.xml','sim_burst');
H2cinjmade6 = readMeta('H2cinj.xml','sim_burst');
L1cinjmade6 = readMeta('L1cinj.xml','sim_burst');

cd /scratch2/xavi/CosmicStrings7

H1ctrigfound7 = readMeta('H1cinjfound.xml','sim_burst');
H2ctrigfound7 = readMeta('H2cinjfound.xml','sim_burst');
L1ctrigfound7 = readMeta('L1cinjfound.xml','sim_burst');

H1cinjmade7 = readMeta('H1cinj.xml','sim_burst');
H2cinjmade7 = readMeta('H2cinj.xml','sim_burst');
L1cinjmade7 = readMeta('L1cinj.xml','sim_burst');

cd /scratch2/xavi/CosmicStrings8

H1ctrigfound8 = readMeta('H1cinjfound.xml','sim_burst');
H2ctrigfound8 = readMeta('H2cinjfound.xml','sim_burst');
L1ctrigfound8 = readMeta('L1cinjfound.xml','sim_burst');

H1cinjmade8 = readMeta('H1cinj.xml','sim_burst');
H2cinjmade8 = readMeta('H2cinj.xml','sim_burst');
L1cinjmade8 = readMeta('L1cinj.xml','sim_burst');


cd /scratch2/xavi/CosmicStrings9

H1ctrigfound9 = readMeta('H1cinjfound.xml','sim_burst');
H2ctrigfound9 = readMeta('H2cinjfound.xml','sim_burst');
L1ctrigfound9 = readMeta('L1cinjfound.xml','sim_burst');

H1cinjmade9 = readMeta('H1cinj.xml','sim_burst');
H2cinjmade9 = readMeta('H2cinj.xml','sim_burst');
L1cinjmade9 = readMeta('L1cinj.xml','sim_burst');

H1foundpeaks=[H1ctrigfound0.hpeak;H1ctrigfound1.hpeak;H1ctrigfound2.hpeak;H1ctrigfound3.hpeak;H1ctrigfound4.hpeak;...
    H1ctrigfound5.hpeak;H1ctrigfound6.hpeak;H1ctrigfound7.hpeak;H1ctrigfound8.hpeak;H1ctrigfound9.hpeak]

H2foundpeaks=[H2ctrigfound0.hpeak;H2ctrigfound1.hpeak;H2ctrigfound2.hpeak;H2ctrigfound3.hpeak;H2ctrigfound4.hpeak;...
    H2ctrigfound5.hpeak;H2ctrigfound6.hpeak;H2ctrigfound7.hpeak;H2ctrigfound8.hpeak;H2ctrigfound9.hpeak]

L1foundpeaks=[L1ctrigfound0.hpeak;L1ctrigfound1.hpeak;L1ctrigfound2.hpeak;L1ctrigfound3.hpeak;L1ctrigfound4.hpeak;...
    L1ctrigfound5.hpeak;L1ctrigfound6.hpeak;L1ctrigfound7.hpeak;L1ctrigfound8.hpeak;L1ctrigfound9.hpeak]


H1madepeaks=[H1cinjmade0.hpeak;H1cinjmade1.hpeak;H1cinjmade2.hpeak;H1cinjmade3.hpeak;H1cinjmade4.hpeak;...
    H1cinjmade5.hpeak;H1cinjmade6.hpeak;H1cinjmade7.hpeak;H1cinjmade8.hpeak;H1cinjmade9.hpeak]

H2madepeaks=[H2cinjmade0.hpeak;H2cinjmade1.hpeak;H2cinjmade2.hpeak;H2cinjmade3.hpeak;H2cinjmade4.hpeak;...
    H2cinjmade5.hpeak;H2cinjmade6.hpeak;H2cinjmade7.hpeak;H2cinjmade8.hpeak;H2cinjmade9.hpeak]

L1madepeaks=[L1cinjmade0.hpeak;L1cinjmade1.hpeak;L1cinjmade2.hpeak;L1cinjmade3.hpeak;L1cinjmade4.hpeak;...
    L1cinjmade5.hpeak;L1cinjmade6.hpeak;L1cinjmade7.hpeak;L1cinjmade8.hpeak;L1cinjmade9.hpeak]

A= -0.2:0.2:2;
nH1found=histc(log10(H1foundpeaks),A);
nH2found=histc(log10(H2foundpeaks),A);
nL1found=histc(log10(L1foundpeaks),A);

nH1made=histc(log10(H1madepeaks),A);
nH2made=histc(log10(H2madepeaks),A);
nL1made=histc(log10(L1madepeaks),A);
 
eH1=nH1found./nH1made;
eH2=nH2found./nH2made;
eL1=nL1found./nL1made;

%figure
plot(A,eH1,'-')
hold on
plot(A,eH1,'bo')
plot(A,eH2,'go')
plot(A,eH2,'g-')
plot(A,eL1,'ro')
plot(A,eL1,'r-')
ylabel('Efficiency')
xlabel('log10(Injected amplitude)')
title('Triple coincident efficiencies for H1 (blue), H2 (green), L1 (red)')

% Find interpolated A50%
for i = 2:size(A,2)
   elH1=eH1(i-1);
   ehH1=eH1(i);
   if elH1 < 0.5
       if ehH1 > 0.5
           A50H1=A(i-1)+(0.5-elH1)*(A(i)-A(i-1))/(ehH1-elH1)
       end
   end

   elH2=eH2(i-1);
   ehH2=eH2(i);
   if elH2 < 0.5
       if ehH2 > 0.5
           A50H2=A(i-1)+(0.5-elH2)*(A(i)-A(i-1))/(ehH2-elH2)
       end
   end
   


   elL1=eL1(i-1);
   ehL1=eL1(i);
   if elL1 < 0.5
       if ehL1 > 0.5
           A50L1=A(i-1)+(0.5-elL1)*(A(i)-A(i-1))/(ehL1-elL1)
       end
   end
  

end

fid=fopen('/home/siemens/lscsoft/lalapps/src/ring/eff_data2.txt','a');
fprintf(fid,'%e %e %e\n',A50H1,A50H2,A50L1);

fclose(fid);

plot(A50H1,0.5,'bx')
plot(A50H2,0.5,'gx')
plot(A50L1,0.5,'rx') 
