clear

t0=10^(17.5);   %age of the universe
xeq=2*10^(7);  %value of x at aeq
c=1;            %Average no of cusps per oscillation
Gamma = 50;     %Gravitational back reaction constant
f=500;           % high frequency cutoff
fl=40;          %low freq cutoff 

p=1;            %Probability of reconnection

Gmu=logspace(-13,-5,1000); %Value of Gmu; this is the variable

n=1;            %Spectral index for small-scale structure

epsilon=1;

alpha=epsilon*(Gamma*Gmu).^n;  %Size of small scale structure

R=1/(365*24*3600);     %Event rate in Hz

b=10^2 * c * alpha.^(-5/3) .* (p*Gamma*Gmu).^(-1) * t0^(-1) * (f*t0)^(-2/3);

x=R./b;

hAlex=(f*t0)^(-1/3) * Gmu.* alpha.^(2/3) ...
     .* x.^(-1/3) .* (1+x).^(-13/33) .* (1+x/xeq).^(3/11);

hMine=6*(fl^(-1/3)-f^(-1/3))*(t0)^(-1/3) * Gmu.* alpha.^(2/3) ...
    .* x.^(-1/3) .* (1+x).^(-13/33) .* (1+x/xeq).^(3/11);


figure(1)
loglog(Gmu,hAlex,'b'); hold on
loglog(Gmu,hMine,'r');
loglog(Gmu,4*1.7*10^(-22),'k:'); hold on
xlabel('G\mu')
ylabel('Strain')

