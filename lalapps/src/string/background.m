clear
path(path,'/home/siemens/ligotools/matlab')

N=50;                %no of slides
tSlide=157080000;  %slide time in nanoseconds

cd /scratch2/xavi/CosmicStrings0/slides/

% This loop finds the minimum and maximum amplitudes
A1=[];
A2=[];
for i=1:N
   
    i
    file_neg=['H1-BURCA_H1L1_M_',num2str(i*tSlide),'_--1-0.xml'];
    file_pos=['H1-BURCA_H1L1_P_',num2str(i*tSlide),'_--1-0.xml'];

    H1ctrig_neg = readMeta(char(file_neg),'sngl_burst');
    H1ctrig_pos = readMeta(char(file_pos),'sngl_burst');
    
    A1=[A1;abs(H1ctrig_neg.amplitude)];
    A2=[A2;abs(H1ctrig_pos.amplitude)];
    
end

Amin=min(min(A1),min(A2))
Amax=max(max(A1),max(A2))

A=Amin:0.5:1.5*Amax;

% This loops finds the histograms
histA1=[];
histA2=[];
for i=1:N
   
    
    file_neg=['H1-BURCA_H1L1_M_',num2str(i*tSlide),'_--1-0.xml'];
    file_pos=['H1-BURCA_H1L1_P_',num2str(i*tSlide),'_--1-0.xml'];

    H1ctrig_neg = readMeta(char(file_neg),'sngl_burst');
    H1ctrig_pos = readMeta(char(file_pos),'sngl_burst');
    
    histA1=[histA1,histc(abs(H1ctrig_neg.amplitude),A)];
    histA2=[histA2,histc(abs(H1ctrig_pos.amplitude),A)];
    
end

histA=[histA1,histA2];

meanhistA=mean(histA,2);
stdhistA=std(histA,0,2); %the 0 is a flag to divide by N-1 

figure
bar(A,meanhistA,'c')
hold on
%errorbar(A,meanhistA,stdhistA,'b.')
plot(A,meanhistA+stdhistA,'bo')
plot(A,meanhistA-stdhistA,'bo')

H1ctrig_foreground=readMeta('H1-BURCA_H1L1_P_0_--1-0.xml','sngl_burst');
histA_foreground=histc(abs(H1ctrig_foreground.amplitude),A);

plot(A,histA_foreground,'rx','MarkerSize',8)
title('Background estimates from slides: Foreground (red x), Background (cyan bars), 1-\sigma (blue o)')
ylabel('Number of occurrences')
xlabel('Amplitude A/10^{-20} [s^{-1/3}]')




