clear
 
t0=10^(17.5);  %age of the universe
zeq=10^(3.94);
aeq=(zeq)^(-4/3) * 2^(-1/3); %value of a at zeq

Ab=9*10^(-21);   %measured amplitude of some event

c=1;           %Average no of cusps per oscillation
Gamma = 50;    %Gravitational back reaction constant

f=75;          %Lowest high frequency cutoff we can detect

p=1;

Gmu=logspace(-13,-5,1000);

n=1;

epsilon=1;

alpha=epsilon*(Gamma*Gmu).^n;

ab=Ab./(t0^(-1/3)*Gmu.*(alpha.^(2/3)));

b=10^2 * c * alpha.^(-5/3) .* (p*Gamma*Gmu).^(-1) * t0^(-1) * (f*t0)^(-2/3);

R=b.* aeq^(33/40) .* ab.^(-11/5) .* ((1+1/aeq*ab).^(33/40)) .* (1+ab).^(-13/8);

theta0=(alpha*f*t0).^(-1/3);
theta=aeq^(3/40)*theta0.* ab.^(-1/5) .* (1+1/aeq*ab).^(3/40) .*  (1+ab).^(1/8);

R=R.*stepfun(1./theta,1);

T=365*24*3600;

figure(3)
loglog(Gmu,R,'b'); hold on
loglog(Gmu,2.303/T,'r')
%loglog(Gmu,1/T,'g')
xlabel('G\mu')
ylabel('Rate[Hz]')
title('One year at LIGO design (A=10^{-21}s^{-1/3}, f_h=150)','FontSize',14)
xlabel('G\mu','FontSize',14)
ylabel('Rate [Hz]','FontSize',14)


