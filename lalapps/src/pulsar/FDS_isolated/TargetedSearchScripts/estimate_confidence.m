close all
clear all

%This script predicts the upper limit of a targeted search given the number
%of SFTs used (per IFO), the weighted noise floor of the SFTs, the largest
%2F value found (and which is used as the cut off in the MC upper limit
%injections, and the declination of the target.

no_H1_sfts = 8770; %Number of H1 sfts
no_H2_sfts = 10005; %Number of H2 sfts
no_L1_sfts = 6782; %Number of L1 sfts

sSh_H1 = 2.56e-22;  %H1 noise level, calculated from SFTs 
sSh_H2 = 4.875e-22; %H2 noise level, calculated from SFTs
sSh_L1 = 4.4e-22; %L1 noise level, calculated from SFTs 
Sh_H1 = sSh_H1^2;
Sh_H2 = sSh_H2^2;
Sh_L1 = sSh_L1^2;

largest_2F_found = 37; %Largest 2F value found in the search

delta = 0.3842252152; %Declination of the Crab pulsar


%%%%%%%%%% End of often changed inputs

cosi = [-1:0.02:1]
psi = [0:10:350]      

%For Hanford - Position information
gamma = 171.8*pi/180;
lambda = 46.45 * pi/180;

%Equations taken from the Appendix of the JKS paper

j3 = (1/128)*(28 - 44*cos(lambda)^2 + 5*sin(2*gamma)^2 * cos(lambda)^4);
j2 = (1/1024)*(68 - 20*cos(lambda)^2 - 13*sin(2*gamma)^2*cos(lambda)^4);
j1 = (1/256)*(4-20*cos(lambda)^2 + 35*sin(2*gamma)^2 * cos(lambda)^4);

e2 = 4*j2-j3*cos(2*delta) + j1*cos(2*delta)^2;
e1 = 4*j1*(cos(delta)^4);

G2 = (1/4)*(1+6*cosi.^2 + cosi.^4);
F2 = (1/4)*(1-cosi.^2).^2;

for x = 1:length(cosi)
	for y = 1:length(psi)
        A2_Han(x,y) = F2(x)*e1*cos(4*psi(y)) + G2(x)*e2;
    end
end

%For Livingston - POsition in formation
gamma = 243.0*pi/180;
lambda = 30.56 * pi/180;

%Equations taken from the Appendix of the JKS paper

j3 = (1/128)*(28 - 44*cos(lambda)^2 + 5*sin(2*gamma)^2 * cos(lambda)^4);
j2 = (1/1024)*(68 - 20*cos(lambda)^2 - 13*sin(2*gamma)^2*cos(lambda)^4);
j1 = (1/256)*(4-20*cos(lambda)^2 + 35*sin(2*gamma)^2 * cos(lambda)^4);

e2 = 4*j2-j3*cos(2*delta) + j1*cos(2*delta)^2;
e1 = 4*j1*(cos(delta)^4);

G2 = (1/4)*(1+6*cosi.^2 + cosi.^4);
F2 = (1/4)*(1-cosi.^2).^2;

for x = 1:length(cosi)
    for y = 1:length(psi)
        A2_Liv(x,y) = F2(x)*e1*cos(4*psi(y)) + G2(x)*e2;
    end
end

%Injection strength values to be calculated
hvalue = [1.2e-24 1.33e-24 1.345e-24 1.36e-24 1.4e-24]

for z = 1:length(hvalue)
    h0 = (hvalue(z)/1.05); %This is a 1.05 factor coming from the exact match vs our real grid spacing with mismatch parameter 0.15.  The upper limits with our real grid spacing were 5% smaller than the upper limits from exact matching.  
    %Effectively, we calculate the expected confidence for a injection strength of 5% smaller than what we actually injected, since the analytical calculations assume perfect matching.

    p = sqrt((A2_Han*h0^2*(no_H2_sfts*1800)/Sh_H2) + (A2_Han*h0^2*(no_H1_sfts*1800)/Sh_H1) + (A2_Liv*h0^2*(no_L1_sfts*1800)/Sh_L1));  %Estimated signal to noise ratio for a given injection strength h0.  Is a matrix covering all polarizations and inclination angles (the A2 parameter).



    for x = 1:101
        for y = 1:36
            q(x,y) = 2*quad(@(F) (0.5*exp(-(2*F+p(x,y)^2)/2).*sqrt((2*F)/p(x,y)^2).*besseli(1,sqrt(2*F*p(x,y)^2))),0,largest_2F_found/2);  %Probability of injection being found as a function of polarization and inclination angle given signal to noise and largest 2F value (which is used as cutoff)
        end
    end

    expected(z) = mean(mean(1-q))  %The expected confidence for a given injection strength, over the entire parameter space.  We have averaged over polarization and inclination angle only at the vert end of the calculation.  
%Effectively, we calculated the probability of being found in 3,636
%discrete boxes in parameter space and assume that when we do a sufficient
%number of injections, this parameter space will be sampled evenly,
%resulting in equal weighting of each box in the final confidence
%calculation.
end

figure(1)
plot( hvalue, expected)
legend('Expected','Location','SouthEast')
xlabel('Injection Strength')
ylabel('Confidence')
title('H1 H2 L1 combined')
