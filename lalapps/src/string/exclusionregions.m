clear

path(path,'/home/siemens/ligotools/matlab')

% Coincident efficiencies

cd /scratch2/xavi/CosmicStringsInjections/

H1ctrigfound0 = readMeta('H1cinjfound.xml','sim_burst');
H2ctrigfound0 = readMeta('H2cinjfound.xml','sim_burst');
L1ctrigfound0 = readMeta('L1cinjfound.xml','sim_burst');

H1cinjmade0 = readMeta('H1cinj.xml','sim_burst');
H2cinjmade0 = readMeta('H2cinj.xml','sim_burst');
L1cinjmade0 = readMeta('L1cinj.xml','sim_burst');

H1foundpeaks=[H1ctrigfound0.hpeak];
H2foundpeaks=[H2ctrigfound0.hpeak];
L1foundpeaks=[L1ctrigfound0.hpeak];

H1madepeaks=[H1cinjmade0.hpeak];
H2madepeaks=[H2cinjmade0.hpeak];
L1madepeaks=[L1cinjmade0.hpeak];

dlnA=0.5;
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

deH1=sqrt(eH1.*(1-eH1)./nH1made);
deH2=sqrt(eH2.*(1-eH2)./nH2made);
deL1=sqrt(eL1.*(1-eL1)./nL1made);

eH1=eH1(1:size(eH1,1)-2);
eH2=eH2(1:size(eH2,1)-2);
eL1=eL1(1:size(eL1,1)-2);

deH1=deH1(1:size(deH1,1)-2);
deH2=deH2(1:size(deH2,1)-2);
deL1=deL1(1:size(deL1,1)-2);

lnA=lnA(1:size(lnA,2)-2);

%Now we can compute the excluded regions as a function of p and n/epsilon.

t0=10^(17.5);  %age of the universe
zeq=10^(3.94);
aeq=(zeq)^(-4/3) * 2^(-1/3); %value of a at zeq

c=1;           %Average no of cusps per oscillation
Gamma = 50;    %Gravitational back reaction constant

f=75;          %Lowest high frequency cutoff we can detect
p=logspace(-4,0,100);        %Probability of reconnection

Gmu=logspace(-8,-5, 4); %Vector with values of Gmu

n=1;%linspace(1,3,100);                    %power law for Gamma*Gmu (as in my PRD paper)
epsilon=logspace(-12,0,100);     %smal-scale structure pre-factor (as in Damour and Vilenkin, 2005 PRD)

A=exp(lnA)*10^(-20);  %turn vector of ln(A) into A; and re-scale by 10^-20 factor

T=1364984;
R90=2.303/T;

figure
loglog(1)
hold on

for i=1:size(Gmu,2)
  for j=1:size(p,2)
    for k=1:size(epsilon,2)      

      alpha=epsilon(k)*(Gamma*Gmu(i))^n;   %value of the size of the small-scale structure
      a=A./(t0^(-1/3)*Gmu(i)*(alpha^(2/3)));   %Reduced amplitude vector      
      b=10^2 * c * alpha^(-5/3)*(p(j)*Gamma*Gmu(i))^(-1) * t0^(-1) * (f*t0)^(-2/3);
      dRdlnA=b*aeq^(33/40)*a.^(-11/5) .* ((1+1/aeq*a).^(33/40)).*(1+a).^(-13/8); %Rate per logarithmic interval of amplitude 

      %apply theta fn cutoff
      theta0=(alpha*f*t0).^(-1/3);
      theta=aeq^(3/40)*theta0.* a.^(-1/5) .* (1+1/aeq*a).^(3/40) .*  (1+a).^(1/8);

      eH1p=eH1.*stepfun(1./theta',1);
      R=dlnA*sum(eH1p.*dRdlnA');
      ul(k,j)=stepfun(R,R90);

    end
  end
  
  contour(p,epsilon,ul,[1 1]);    
  xlabel('Reconnection  probability p','FontSize',14)
  ylabel('\epsilon (small-scale structure parameter)','FontSize',14)
  title(['Exclusion regions',num2str(Gmu(i))],'FontSize',14)

end
