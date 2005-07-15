clear

path(path,'/home/siemens/ligotools/matlab')
%cd /scratch3/xavi/CosmicStringInjections

H1trigfound = readMeta('H1injfound.xml','sim_burst');
H2trigfound = readMeta('H2injfound.xml','sim_burst');
L1trigfound = readMeta('L1injfound.xml','sim_burst');

H1injmade = readMeta('H1inj.xml','sim_burst');
H2injmade = readMeta('H2inj.xml','sim_burst');
L1injmade = readMeta('L1inj.xml','sim_burst');

x= 0.1:0.2:2;
nH1found=histc(log10(H1trigfound.hpeak),x);
nH2found=histc(log10(H2trigfound.hpeak),x);
nL1found=histc(log10(L1trigfound.hpeak),x);

nH1made=histc(log10(H1injmade.hpeak),x);
nH2made=histc(log10(H2injmade.hpeak),x);
nL1made=histc(log10(L1injmade.hpeak),x);
 
eH1=nH1found./nH1made;
eH2=nH2found./nH2made;
eL1=nL1found./nL1made;

figure
plot(x,eH1,'o')
hold on
grid on
plot(x,eH1,'-')
plot(x,eH2,'go')
plot(x,eH2,'g-')
plot(x,eL1,'ro')
plot(x,eL1,'r-')
ylabel('Efficiency')
xlabel('log10(Injected amplitude)')
title('Single IFO efficiencies for H1 (blue), H2 (green), L1 (red); triple coincident (dashed)')

% Coincident efficiencies
 
H1ctrigfound = readMeta('H1cinjfound.xml','sim_burst');
H2ctrigfound = readMeta('H2cinjfound.xml','sim_burst');
L1ctrigfound = readMeta('L1cinjfound.xml','sim_burst');

H1cinjmade = readMeta('H1cinj.xml','sim_burst');
H2cinjmade = readMeta('H2cinj.xml','sim_burst');
L1cinjmade = readMeta('L1cinj.xml','sim_burst');

nH1found=histc(log10(H1ctrigfound.hpeak),x);
nH2found=histc(log10(H2ctrigfound.hpeak),x);
nL1found=histc(log10(L1ctrigfound.hpeak),x);

nH1made=histc(log10(H1cinjmade.hpeak),x);
nH2made=histc(log10(H2cinjmade.hpeak),x);
nL1made=histc(log10(L1cinjmade.hpeak),x);
 
eH1=nH1found./nH1made;
eH2=nH2found./nH2made;
eL1=nL1found./nL1made;

plot(x,eH1,'--')
plot(x,eH2,'go')
plot(x,eH2,'g--')
plot(x,eL1,'ro')
plot(x,eL1,'r--')


 
