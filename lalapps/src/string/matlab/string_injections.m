clear
path(path,'/home/siemens/ligotools/matlab')

injfile1='HL-INJECTIONS-795169195-1804.xml';
injfile2='HL-INJECTIONS-795171034-685.xml';
injfile3='HL-INJECTIONS-795171754-2710.xml';
injfile4='HL-INJECTIONS-795174556-3274.xml';
injfile5='HL-INJECTIONS-795177865-5056.xml';
injfile6='HL-INJECTIONS-795183281-4288.xml';
injfile7='HL-INJECTIONS-795187604-7381.xml';
injfile8='HL-INJECTIONS-795195020-5516.xml';
injfile9='HL-INJECTIONS-795200571-8206.xml';
injfile10='HL-INJECTIONS-795208812-2792.xml';

H1trigfile='H1triggers.xml';
H2trigfile='H2triggers.xml';
L1trigfile='L1triggers.xml';
H1ctrigfile='H1ctriggers.xml';
H2ctrigfile='H2ctriggers.xml';
L1ctrigfile='L1ctriggers.xml';

cd /home/siemens/lscsoft/lalapps/src/ring/CShtInjectionsLONG

inj1=readMeta(injfile1,'sim_burst');
inj2=readMeta(injfile2,'sim_burst');
inj3=readMeta(injfile3,'sim_burst');
inj4=readMeta(injfile4,'sim_burst');
inj5=readMeta(injfile5,'sim_burst');
inj6=readMeta(injfile6,'sim_burst');
inj7=readMeta(injfile7,'sim_burst');
inj8=readMeta(injfile8,'sim_burst');
inj9=readMeta(injfile9,'sim_burst');
inj10=readMeta(injfile10,'sim_burst');

H1trig = readMeta(H1trigfile,'sngl_burst');
H2trig = readMeta(H2trigfile,'sngl_burst');
L1trig = readMeta(L1trigfile,'sngl_burst');

H1ctrig = readMeta(H1ctrigfile,'sngl_burst');
H2ctrig = readMeta(H2ctrigfile,'sngl_burst');
L1ctrig = readMeta(L1ctrigfile,'sngl_burst');

H1trig_inj = readMeta('H1out.xml','sngl_burst');
H2trig_inj = readMeta('H2out.xml','sngl_burst');
L1trig_inj = readMeta('L1out.xml','sngl_burst');

H1injfound = readMeta('H1injfound.xml','sim_burst');
H2injfound = readMeta('H2injfound.xml','sim_burst');
L1injfound = readMeta('L1injfound.xml','sim_burst');

H1ctrig_inj = readMeta('H1cout.xml','sngl_burst');
H2ctrig_inj = readMeta('H2cout.xml','sngl_burst');
L1ctrig_inj = readMeta('L1cout.xml','sngl_burst');

H1cinjfound = readMeta('H1cinjfound.xml','sim_burst');
H2cinjfound = readMeta('H2cinjfound.xml','sim_burst');
L1cinjfound = readMeta('L1cinjfound.xml','sim_burst');

gps = 795169179;

injgeocent_peak_time=[inj1.geocent_peak_time;inj2.geocent_peak_time;inj3.geocent_peak_time;inj4.geocent_peak_time;inj5.geocent_peak_time;inj6.geocent_peak_time; inj7.geocent_peak_time;inj8.geocent_peak_time;inj9.geocent_peak_time;inj10.geocent_peak_time];

injgeocent_peak_time_ns=[inj1.geocent_peak_time_ns;inj2.geocent_peak_time_ns;inj3.geocent_peak_time_ns; inj4.geocent_peak_time_ns;inj5.geocent_peak_time_ns;inj6.geocent_peak_time_ns; inj7.geocent_peak_time_ns;inj8.geocent_peak_time_ns;inj9.geocent_peak_time_ns; inj10.geocent_peak_time_ns];
  
injhpeak=[inj1.hpeak;inj2.hpeak;inj3.hpeak;inj4.hpeak;inj5.hpeak;inj6.hpeak;inj7.hpeak;inj8.hpeak;inj9.hpeak;inj10.hpeak];

figure
semilogy((injgeocent_peak_time+injgeocent_peak_time_ns*1e-9)-gps,injhpeak,'rx'); grid on
hold on;
plot(H1trig.peak_time+H1trig.peak_time_ns*1e-9-gps,abs(H1trig.amplitude),'bo','MarkerSize',8)
plot(H2trig.peak_time+H2trig.peak_time_ns*1e-9-gps,abs(H2trig.amplitude),'gd','MarkerSize',8)
plot(L1trig.peak_time+L1trig.peak_time_ns*1e-9-gps,abs(L1trig.amplitude),'rd','MarkerSize',8)
plot((injgeocent_peak_time+injgeocent_peak_time_ns*1e-9)-gps,injhpeak,'k.','MarkerSize',10)
title('Trigger amplitudes: H1 (blue), H1 (green), L1 (red), injected (black)','FontSize',12)
xlabel('Time in seconds','FontSize',12)
ylabel('Amplitude (10^{-20} s^{-1/3})','FontSize',12)

figure
semilogy((injgeocent_peak_time+injgeocent_peak_time_ns*1e-9)-gps,injhpeak,'kx'); grid on
hold on;
plot(H1injfound.geocent_peak_time+H1injfound.geocent_peak_time_ns*1e-9-gps,abs(H1injfound.hpeak),'bo','MarkerSize',8)
plot(H2injfound.geocent_peak_time+H2injfound.geocent_peak_time_ns*1e-9-gps,abs(H2injfound.hpeak),'gd','MarkerSize',8)
plot(L1injfound.geocent_peak_time+L1injfound.geocent_peak_time_ns*1e-9-gps,abs(L1injfound.hpeak),'rd','MarkerSize',8)
semilogy((injgeocent_peak_time+injgeocent_peak_time_ns*1e-9)-gps,injhpeak,'k.'); grid on
title('Single IFO injections found: H1 (blue), H1 (green), L1 (red), injected (black)','FontSize',12)
xlabel('Time in seconds','FontSize',12)
ylabel('Injected amplitude (10^{-20} s^{-1/3})','FontSize',12)


figure
semilogy((injgeocent_peak_time+injgeocent_peak_time_ns*1e-9)-gps,injhpeak,'rx'); grid on
hold on;
plot(H1ctrig.peak_time+H1ctrig.peak_time_ns*1e-9-gps,abs(H1ctrig.amplitude),'bo','MarkerSize',8)
plot(H2ctrig.peak_time+H2ctrig.peak_time_ns*1e-9-gps,abs(H2ctrig.amplitude),'gd','MarkerSize',8)
plot(L1ctrig.peak_time+L1ctrig.peak_time_ns*1e-9-gps,abs(L1ctrig.amplitude),'rd','MarkerSize',8)
plot((injgeocent_peak_time+injgeocent_peak_time_ns*1e-9)-gps,injhpeak,'k.','MarkerSize',10)
title('Post coincidence trigger amplitudes: H1 (blue), H1 (green), L1 (red), injected (black)','FontSize',12)
xlabel('Time in seconds','FontSize',12)
ylabel('Amplitude (10^{-20} s^{-1/3})','FontSize',12)

figure
semilogy((injgeocent_peak_time+injgeocent_peak_time_ns*1e-9)-gps,injhpeak,'kx'); grid on
hold on;
plot(H1cinjfound.geocent_peak_time+H1cinjfound.geocent_peak_time_ns*1e-9-gps,abs(H1cinjfound.hpeak),'bo','MarkerSize',8)
plot(H2cinjfound.geocent_peak_time+H2cinjfound.geocent_peak_time_ns*1e-9-gps,abs(H2cinjfound.hpeak),'gd','MarkerSize',8)
plot(L1cinjfound.geocent_peak_time+L1cinjfound.geocent_peak_time_ns*1e-9-gps,abs(L1cinjfound.hpeak),'rd','MarkerSize',8)
semilogy((injgeocent_peak_time+injgeocent_peak_time_ns*1e-9)-gps,injhpeak,'k.'); grid on
title('Injections found in coincidence: H1 (blue), H1 (green), L1 (red), injected (black)','FontSize',12)
xlabel('Time in seconds','FontSize',12)
ylabel('Injected amplitude (10^{-20} s^{-1/3})','FontSize',12)

a=(H1injfound.hpeak-(H1trig_inj.amplitude))./H1injfound.hpeak;
b=(H2injfound.hpeak-(H2trig_inj.amplitude))./H2injfound.hpeak;
c=(L1injfound.hpeak-(L1trig_inj.amplitude))./L1injfound.hpeak;

figure
loglog(H1injfound.hpeak,abs(a),'bo','MarkerSize',8); grid on
hold on
plot(H2injfound.hpeak,abs(b),'go','MarkerSize',8); grid on
plot(L1injfound.hpeak,abs(c),'ro','MarkerSize',8); grid on
title('Relative errors in recovered amplitudes: H1 (blue), H1 (green), L1 (red)','FontSize',12)
xlabel('Injected amplitudes(10^{-20} s^{-1/3})','FontSize',12)
ylabel('Relative difference (A_i-A_d)/A_i','FontSize',12)

figure
loglog(1./H1trig_inj.snr,abs(a),'bo','MarkerSize',8); grid on
hold on
plot(1./H2trig_inj.snr,abs(b),'go','MarkerSize',8); grid on
plot(1./L1trig_inj.snr,abs(c),'ro','MarkerSize',8); grid on
x=0.01:0.01:1;
y=x;
plot(x,y,'m','LineWidth',2)
title('Errors in amplitudes and SNR: H1 (blue), H1 (green), L1 (red)','FontSize',12)
xlabel('1/SNR','FontSize',12)
ylabel('Relative difference (A_i-A_d)/A_i','FontSize',12)


x= -5:0.01:5;
na=histc(a,x);
nb=histc(b,x);
nc=histc(c,x);
figure
bar(x,na,'b')
title('Histogram of relative errors in the amplitude: H1 (blue), H2 (green), L1 (red)','FontSize',12)
xlabel('relative difference (Ai-Ad)/Ai','FontSize',12)
ylabel('Number of occurences','FontSize',12)
hold on
bar(x,nb,'g')
bar(x,nc,'r')

at=(H1injfound.geocent_peak_time+H1injfound.geocent_peak_time_ns*1e-9)-(H1trig_inj.peak_time+H1trig_inj.peak_time_ns*1e-9);
bt=(H2injfound.geocent_peak_time+H2injfound.geocent_peak_time_ns*1e-9)-(H2trig_inj.peak_time+H2trig_inj.peak_time_ns*1e-9);
ct=(L1injfound.geocent_peak_time+L1injfound.geocent_peak_time_ns*1e-9)-(L1trig_inj.peak_time+L1trig_inj.peak_time_ns*1e-9);

x= -1e-2:1e-5:1e-2;
nat=histc(at,x);
nbt=histc(bt,x);
nct=histc(ct,x);
figure
bar(x,nat,'b')
title('Histogram of injected - detected peak time H1 (blue),  H2 (green), L1 (red)','FontSize',12)
xlabel('relative difference ti-td in seconds','FontSize',12)
ylabel('Number of occurences','FontSize',12)
hold on
bar(x,nbt,'g')
bar(x,nct,'r')

figure
plot(at,a,'bo')
hold on
plot(bt,b,'gd')
plot(ct,c,'rd')
title('Errors in amplitudes vs differences in trigger peak time: H1 (blue),  H2 (green), L1 (red)','FontSize',12)
xlabel('Difference ti-td in seconds','FontSize',12)
ylabel('Relative difference in amplitude','FontSize',12)





 
 
 
