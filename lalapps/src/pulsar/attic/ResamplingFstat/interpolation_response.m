% Script to compute the magntude response of an interpolating function

band=5.0; % replace B with the band in Hz of your heterodyned time series
df=1/10000;  % replace DF with the frequency resolution you want 
        % (best is to use frequency resolution of F-stat i.e. 1/T)

f=0:df:band;

b=[0.5 0.5]

% This is where the magic happens. The numbers b1 b2 etc are the 
% coefficients used in the interpolation. For example, for linear 
% interpolation they would be b1=0.5, b2=0.5, i.e.:
% b=[0.5 0.5];

a=[1 zeros(size(b,2)-1,1)'];

interp_filter=tf(b,a,1/(2*band));
hI=squeeze(freqresp(interp_filter,2*pi*f));

%figure
plot(f,abs(hI),'b.')
xlabel('Hz','FontSize',12)
ylabel('Magnitude','FontSize',12)
title('Magnitude response','FontSize',12)
