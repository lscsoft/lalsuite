clear

path(path,'/home/siemens/ligotools/matlab')

% Coincident efficiencies

cd /scratch2/xavi/CosmicStrings0/

H1ctrigfound0 = readMeta('H1cinjfound.xml','sim_burst');
H2ctrigfound0 = readMeta('H2cinjfound.xml','sim_burst');
L1ctrigfound0 = readMeta('L1cinjfound.xml','sim_burst');

H1cinjmade0 = readMeta('H1cinj.xml','sim_burst');
H2cinjmade0 = readMeta('H2cinj.xml','sim_burst');
L1cinjmade0 = readMeta('L1cinj.xml','sim_burst');

% gps = 793154935;
% figure
% semilogy((H1cinjmade0.geocent_peak_time+H1cinjmade0.geocent_peak_time_ns*1e-9)-gps,H1cinjmade0.hpeak,'rx'); grid on
% hold on;
% semilogy((H1ctrigfound0.geocent_peak_time+H1ctrigfound0.geocent_peak_time_ns*1e-9)-gps,H1ctrigfound0.hpeak,'ko'); grid on

H1foundpeaks=[H1ctrigfound0.hpeak];
H2foundpeaks=[H2ctrigfound0.hpeak];
L1foundpeaks=[L1ctrigfound0.hpeak];

H1madepeaks=[H1cinjmade0.hpeak];
H2madepeaks=[H2cinjmade0.hpeak];
L1madepeaks=[L1cinjmade0.hpeak];

dlnA=0.7;
lnA=-0.46:dlnA:7;

nH1found=histc(log(H1foundpeaks),lnA);
nH2found=histc(log(H2foundpeaks),lnA);
nL1found=histc(log(L1foundpeaks),lnA);

nH1made=histc(log(H1madepeaks),lnA);
nH2made=histc(log(H2madepeaks),lnA);
nL1made=histc(log(L1madepeaks),lnA);
 
eH1=nH1found./nH1made;
eH2=nH2found./nH2made;
eL1=nL1found./nL1made;

eH1=eH1(1:size(eH1,1)-2);
eH2=eH2(1:size(eH2,1)-2);
eL1=eL1(1:size(eL1,1)-2);

lnA=lnA(1:size(lnA,2)-2);

figure
plot(lnA,eH1,'-')
hold on
plot(lnA,eH1,'bo')
plot(lnA,eH2,'go')
plot(lnA,eH2,'g-')
plot(lnA,eL1,'ro')
plot(lnA,eL1,'r-')
plot(lnA,eH1,'-')
plot(lnA,eH1,'bo')
ylabel('\epsilon (A^*,A)','FontSize',12)
xlabel('ln(A)','FontSize',12)
title('Triple coincident efficiencies for H1 (blue), H2 (green), L1 (red)','FontSize',12)


%Now we can compute the excluded regions as a function of Gmu

t0=10^(17.5);  %age of the universe
zeq=10^(3.94);
aeq=(zeq)^(-4/3) * 2^(-1/3); %value of a at zeq

c=1;           %Average no of cusps per oscillation
Gamma = 50;    %Gravitational back reaction constant

f=75;          %Lowest high frequency cutoff we can detect
p=1e-6;           %Probability of reconnection

Gmu=logspace(-13,-5,8); %Vector with values of Gmu

n=1;           %power law for Gamma*Gmu (as in my PRD paper)
epsilon=1;     %smal-scale structure pre-factor (as in Damour and Vilenkin, 2005 PRD)

A=exp(lnA)*10^(-20);  %turn vector of ln(A) back into A; and re-scale by 10^-20 factor

for i=1:size(Gmu,2)
        
    alpha=epsilon*(Gamma*Gmu(i))^n;   %value of the size of teh small-scale structure
    a=A./(t0^(-1/3)*Gmu(i)*(alpha^(2/3)));   %Reduced amplitude vector      
    b=10^2 * c * alpha^(-5/3)*(p*Gamma*Gmu(i))^(-1) * t0^(-1) * (f*t0)^(-2/3);
    dRdlnA=b*aeq^(33/40)*a.^(-11/5) .* ((1+1/aeq*a).^(33/40)).*(1+a).^(-13/8); %Rate per logarithmic interval of amplitude 
    R(i)=dlnA*sum(eH1.*dRdlnA');

    R(i)=R(i)+dRdlnA(size(dRdlnA,2));   %This last addition is for events with amplitude larger than the largest amplitude 
                             %(which I assume we detect with efficiency 1)
end


T=204211;

figure(2)
loglog(Gmu,R,'b'); hold on
loglog(Gmu,2.303/T,'r')
xlabel('G\mu')
ylabel('Rate[Hz]')
title('Playground upper limit','FontSize',14)
xlabel('G\mu','FontSize',14)
ylabel('Rate [Hz]','FontSize',14)
