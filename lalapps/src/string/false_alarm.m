% clear

path(path,'/home/siemens/ligotools/matlab')

cd /scratch2/xavi/CosmicStrings

H1trigfile='H1triggers.xml'
H2trigfile='H2triggers.xml'
L1trigfile='L1triggers.xml'
H1ctrigfile='H1ctriggers.xml'
H2ctrigfile='H2ctriggers.xml'
L1ctrigfile='L1ctriggers.xml'
% 
% H1trig = readMeta(H1trigfile,'sngl_burst');
% H2trig = readMeta(H2trigfile,'sngl_burst');
% L1trig = readMeta(L1trigfile,'sngl_burst');
H1ctrig = readMeta(H1ctrigfile,'sngl_burst');
H2ctrig = readMeta(H2ctrigfile,'sngl_burst');
L1ctrig = readMeta(L1ctrigfile,'sngl_burst');

T=41872;

fid=fopen('/home/siemens/lscsoft/lalapps/src/ring/fa_data.txt','w');

for i=1:400

    thres=3.75+(i-1)*0.25;
    %H1:
    [j,k]=find(H1trig.snr>=thres);
    RH1=size(H1trig.snr(j),1)/T;
    [j,k]=find(H1ctrig.snr>=thres);
    RcH1=size(H1ctrig.snr(j),1)/T;

    %L1:
    [j,k]=find(L1trig.snr>=thres);
    RL1=size(L1trig.snr(j),1)/T;
    [j,k]=find(L1ctrig.snr>=thres);
    RcL1=size(L1ctrig.snr(j),1)/T;

    %H2:
    [j,k]=find(H2trig.snr>=thres);
    RH2=size(H2trig.snr(j),1)/T;
    [j,k]=find(H2ctrig.snr>=thres);
    RcH2=size(H2ctrig.snr(j),1)/T;

    FA(i,1)=thres;FA(i,2)=RL1;FA(i,3)=RH1;FA(i,4)=RH2;FA(i,5)=RcL1;FA(i,6)=RcH1;FA(i,7)=RcH2;
    FA(i,8)=FA(i,2)*FA(i,3)*FA(i,4)*(2*0.004)*(2*0.020);
    FA(i,9)=FA(i,2)*FA(1,3)*FA(1,4)*(2*0.004)*(2*0.020);
    FA(i,10)=FA(1,2)*FA(i,3)*FA(1,4)*(2*0.004)*(2*0.020);
    FA(i,11)=FA(1,2)*FA(1,3)*FA(i,4)*(2*0.004)*(2*0.020);
    fprintf(fid,'%e %e %e %e %e %e %e %e %e %e %e\n',FA(i,1),FA(i,2),FA(i,3),FA(i,4),FA(i,5),FA(i,6),FA(i,7),FA(i,8),FA(i,9),FA(i,10),FA(i,11));
end
fclose(fid);
figure
semilogy(FA(:,1),FA(:,2),'r-')
hold on
semilogy(FA(:,1),FA(:,2),'ro')
semilogy(FA(:,1),FA(:,3),'bo')
semilogy(FA(:,1),FA(:,3),'b-')
semilogy(FA(:,1),FA(:,4),'g-')
semilogy(FA(:,1),FA(:,4),'go')
semilogy(FA(:,1),FA(:,5),'ro')
semilogy(FA(:,1),FA(:,5),'r-')
semilogy(FA(:,1),FA(:,6),'b-')
semilogy(FA(:,1),FA(:,6),'bo')
semilogy(FA(:,1),FA(:,7),'go')
semilogy(FA(:,1),FA(:,7),'g-')
semilogy(FA(:,1),FA(:,9),'r--')
semilogy(FA(:,1),FA(:,10),'b--')
semilogy(FA(:,1),FA(:,11),'g--')
semilogy(FA(:,1),FA(:,8),'k:')
ylabel('Event rate (Hz)')
xlabel('SNR')
title('Single IFO and triple coincidence false alarm rates for H1 (blue), H2 (green), L1 (red)')


tH1=H1trig.peak_time+1e-9*H1trig.peak_time_ns;
tH2=H2trig.peak_time+1e-9*H2trig.peak_time_ns;
tL1=L1trig.peak_time+1e-9*L1trig.peak_time_ns;

tH1=sort(tH1);
tH2=sort(tH2);
tL1=sort(tL1);

x=0:0.01:10;
yH1=histc(diff(tH1),x)/0.01/size(diff(tH1),1);
yH2=histc(diff(tH2),x)/0.01/size(diff(tH2),1);
yL1=histc(diff(tL1),x)/0.01/size(diff(tL1),1);

figure
semilogy(x,yH1,'b')
hold on
semilogy(x,yH2,'g')
semilogy(x,yL1,'r')
y=FA(1,3)*exp(-FA(1,3)*x);
plot(x,y)

