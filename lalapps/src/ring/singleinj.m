
function singleinj(found_list, ifo)

% 
% NULL = singleinj ( found_list, ifo )
%
% Generate parameter estimation plots from found_list, which is a column vector of 
% filenames output from lalapps_coincringread.  ifo is a string of the name of the
% interferometer for which the found parameters will be compared to the injected. 
% 
% plots are saved in injdetplots/ which is generated if it does not already exist.
%
% EXAMPLE:
% 
% found_list = ['injH1H2L1coincs_2001.xml';'injH1H2L1coincs_2002.xml'];
% singleinj(found_list, 'H1')
%
% Lisa Goggin, Mar 15th 2006
% modified P. T. Baker 30 June 2009

  % read in the injection lists
  injH1t=load('injH1_trip.mat');
  injL1t=load('injL1_trip.mat');
  injH2t=load('injH2_trip.mat');
  injH1d=load('injH1d.mat');
  injH2d=load('injH2d.mat');
  injL1d=load('injL1d.mat');
  injH1inL1d=load('injH1inL1_doub.mat');
  injL1inH1d=load('injL1inH1_doub.mat');
  injH1inH2d=load('injH1inH2_doub.mat');
  injH2inH1d=load('injH2inH1_doub.mat');
  injL1inH2d=load('injL1inH2_doub.mat');
  injH2inL1d=load('injH2inL1_doub.mat');

  % if no cut
  injH1inL1dcut.ind = [];
  injH1inH2dcut.ind = [];
  injH2inL1dcut.ind = [];
  injH1inL1daftercut.ind = injH1inL1d.ind;
  injH1inH2daftercut.ind = injH1inH2d.ind;
  injH2inL1daftercut.ind = injH2inL1d.ind;


  %=========================== Injected Quantities ===========================%

  % Read found sim_ringdown tables
  N_files = length(found_list(:,1))

  eval(['fsim=readMeta(found_list(1,:),''sim_ringdown'',0,''h_start_time,h_start_time_ns,l_start_time,l_start_time_ns,mass,spin,frequency,quality,eff_dist_h,eff_dist_l,distance,hrss_h,hrss_l'');'])

  A=1:N_files;
  for k=1:length(fsim.mass)
    fsim.run(k)=A(1);
  end
  fsim.run=transpose(fsim.run);

  for i=2:N_files
    eval(['fsimi=readMeta(found_list(i,:),''sim_ringdown'',0,''h_start_time,h_start_time_ns,l_start_time,l_start_time_ns,mass,spin,quality,frequency,eff_dist_h,eff_dist_l,distance,hrss_h,hrss_l'');'])

    for k=1:length(fsimi.quality)
      fsimi.run(k)=A(i);
    end
    fsimi.run=transpose(fsimi.run);

    fsim.h_start_time=[fsim.h_start_time;fsimi.h_start_time];
    fsim.h_start_time_ns=[fsim.h_start_time_ns;fsimi.h_start_time_ns];
    fsim.l_start_time=[fsim.l_start_time;fsimi.l_start_time];
    fsim.l_start_time_ns=[fsim.l_start_time_ns;fsimi.l_start_time_ns];
    fsim.frequency=[fsim.frequency;fsimi.frequency];
    fsim.quality=[fsim.quality;fsimi.quality];
    fsim.mass=[fsim.mass;fsimi.mass];  
    fsim.spin=[fsim.spin;fsimi.spin];
    fsim.eff_dist_h=[fsim.eff_dist_h;fsimi.eff_dist_h];
    fsim.eff_dist_l=[fsim.eff_dist_l;fsimi.eff_dist_l];
    fsim.distance=[fsim.distance;fsimi.distance];
    fsim.run=[fsim.run;fsimi.run];
    fsim.hrss_h=[fsim.hrss_h;fsimi.hrss_h];
    fsim.hrss_l=[fsim.hrss_l;fsimi.hrss_l];

  end

  fsim.th=fsim.h_start_time+fsim.h_start_time_ns/1e9;
  fsim.tl=fsim.l_start_time+fsim.l_start_time_ns/1e9;
  fsim.f=fsim.frequency;
  fsim.q=fsim.quality;
  fsim.m=fsim.mass;
  fsim.a=fsim.spin;
  fsim.d=fsim.distance;
  fsim.dh=fsim.eff_dist_h;
  fsim.dl=fsim.eff_dist_l;

  for i=1:length(fsim.th)
    fsim.ind(i)=i;
  end

  % the index of fsim is the same as the fsim.ind
  fsimH1L1d.ind=fsim.ind(injH1inL1daftercut.ind);
  fsimH1L1d.f=fsim.f(injH1inL1daftercut.ind);
  fsimH1L1d.d=fsim.d(injH1inL1daftercut.ind);
  fsimH1L1d.dh=fsim.dh(injH1inL1daftercut.ind);
  fsimH1L1d.dl=fsim.dl(injH1inL1daftercut.ind);
  fsimH1L1d.th=fsim.th(injH1inL1daftercut.ind);
  fsimH1L1d.tl=fsim.tl(injH1inL1daftercut.ind);
  fsimH1L1d.run=fsim.run(injH1inL1daftercut.ind);
  fsimH1L1d.hrss_h=fsim.hrss_h(injH1inL1daftercut.ind);
  fsimH1L1d.hrss_l=fsim.hrss_l(injH1inL1daftercut.ind);

  fsimH1H2d.ind=fsim.ind(injH1inH2daftercut.ind);
  fsimH1H2d.f=fsim.f(injH1inH2daftercut.ind);
  fsimH1H2d.d=fsim.d(injH1inH2daftercut.ind);
  fsimH1H2d.dh=fsim.dh(injH1inH2daftercut.ind);
  fsimH1H2d.dl=fsim.dl(injH1inH2daftercut.ind);
  fsimH1H2d.th=fsim.th(injH1inH2daftercut.ind);
  fsimH1H2d.tl=fsim.tl(injH1inH2daftercut.ind);
  fsimH1H2d.run=fsim.run(injH1inH2daftercut.ind);
  fsimH1H2d.hrss_h=fsim.hrss_h(injH1inH2daftercut.ind);
  fsimH1H2d.hrss_l=fsim.hrss_l(injH1inH2daftercut.ind);

  fsimL1H2d.ind=fsim.ind(injH2inL1daftercut.ind);
  fsimL1H2d.f=fsim.f(injH2inL1daftercut.ind);
  fsimL1H2d.d=fsim.d(injH2inL1daftercut.ind);
  fsimL1H2d.dh=fsim.dh(injH2inL1daftercut.ind);
  fsimL1H2d.dl=fsim.dl(injH2inL1daftercut.ind);
  fsimL1H2d.th=fsim.th(injH2inL1daftercut.ind);
  fsimL1H2d.tl=fsim.tl(injH2inL1daftercut.ind);
  fsimL1H2d.run=fsim.run(injH2inL1daftercut.ind);
  fsimL1H2d.hrss_h=fsim.hrss_h(injH2inL1daftercut.ind);
  fsimL1H2d.hrss_l=fsim.hrss_l(injH2inL1daftercut.ind);

  foundd=[fsimH1L1d.ind,fsimH1H2d.ind,fsimL1H2d.ind];
  threshcut.ind=fsim.ind([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);
  threshcut.f=fsim.f([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);
  threshcut.d=fsim.d([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);
  threshcut.dh=fsim.dh([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);
  threshcut.dl=fsim.dl([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);
  threshcut.th=fsim.th([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);
  threshcut.tl=fsim.tl([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);
  threshcut.run=fsim.run([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);
  threshcut.hrss_h=fsim.hrss_h([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);
  threshcut.hrss_l=fsim.hrss_l([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);

  % get site for input ifo
  if strcmp(ifo,'H1')
    site='h';
  elseif strcmp(ifo,'H2')
    site='h';
  elseif strcmp(ifo,'L1')
    site='l';
  end

  % decay time tau=Q/(pi*f)
  fsim.tau=fsim.q./pi./fsim.f;
  eval(['inj' ifo 't.tau=inj' ifo 't.Q./pi./inj' ifo 't.f;'])
  
  % transpose fsim vectors
  fsim.th=fsim.th';
  fsim.tl=fsim.tl';
  fsim.dh=fsim.dh';
  fsim.dl=fsim.dl';
  fsim.f=fsim.f';
  fsim.q=fsim.q';
  fsim.a=fsim.a';
  fsim.m=fsim.m';
  fsim.d=fsim.d';

  if ~exist('injdetplots','dir')
    mkdir injdetplots
  end
  cd injdetplots

  %%%%%%%%%%%%%%%%%%%%%%% difference scatter plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%% TIME %%%%%%%%%%%%%%%
  % t_det-t_inj v effective distance
  figure
  eval(['semilogx(fsim.d' site '(inj' ifo 't.ind),inj' ifo 't.t-fsim.t' site '(inj' ifo 't.ind),''.'')'])
  grid on
  xlabel('d_{eff_{inj}} / Mpc')
  ylabel('t_{det}-t_{inj} / s')
  eval(['title(''' ifo ': Detected-injected time versus injected effective distance'')'])
  eval(['saveas(gcf,''' ifo 'injdettvd.png'')'])

  % t_det-t_inj v tau
  figure
  eval(['semilogx(fsim.tau(inj' ifo 't.ind),inj' ifo 't.t-fsim.t' site '(inj' ifo 't.ind),''.'')'])
  grid on
  xlabel('\tau / s')
  ylabel('t_{det}-t_{inj} / s')
  eval(['title(''' ifo ': Detected-injected time versus injected \tau=Q/pi/f'')'])
  eval(['saveas(gcf,''' ifo 'injdettvtau.png'')'])

  % t_det-t_inj v Q
  figure
  eval(['plot(fsim.q(inj' ifo 't.ind),inj' ifo 't.t-fsim.t' site '(inj' ifo 't.ind),''.'')'])
  grid on
  xlabel('Q')
  ylabel('t_{det}-t_{inj} / s')
  eval(['title(''' ifo ': Detected-injected time versus injected quality'')'])
  eval(['saveas(gcf,''' ifo 'injdettvQ.png'')'])

  % (t_det-t_inj)*1000 v f 
  figure
  h_p=eval(['semilogx(fsim.f(inj' ifo 't.ind),1000*(inj' ifo 't.t-fsim.t' site '(inj' ifo 't.ind)),''.'')']);
  grid on
  h_xlab=xlabel('f / Hz');
  h_ylab=ylabel('t_{det}-t_{inj} / ms');
  grid on
  set(h_xlab,'FontSize',16,'FontName','Times');
  set(h_ylab,'FontSize',16,'FontName','Times');
  set(gca,'FontSize',16,'FontName','Times');
  set(h_p,'MarkerSize',6)
  %eval(['title(''' ifo ': Detected-injected time versus injected frequency'')'])
  eval(['saveas(gcf,''' ifo 'injdettvf.png'')'])

  % histogram of det-inj time
  figure
  eval(['hist(inj' ifo 't.t-fsim.t' site '(inj' ifo 't.ind),100)'])
  grid on
  xlabel('t_{det}-t_{inj} / s')
  eval(['title(''' ifo ': Detected-injected time'')'])
  eval(['saveas(gcf,''' ifo 'histinjdett.png'')'])

  % all three ifos
  [aL1,bL1]=hist((injL1t.t-fsim.tl(injL1t.ind))*1000,100);
  [aH1,bH1]=hist((injH1t.t-fsim.th(injH1t.ind))*1000,bL1);
  [aH2,bH2]=hist((injH2t.t-fsim.th(injH2t.ind))*1000,bL1);

  figure
  h_H1=plot(bH1,aH1,'r');
  hold on
  h_H2=plot(bH2,aH2,'b');
  h_L1=plot(bL1,aL1,'g');
  set(h_H1,'LineWidth',2)
  set(h_H2,'LineWidth',2)
  set(h_L1,'LineWidth',2)
  h_leg=legend('H1','H2','L1');
  h_xlab=xlabel('\deltat / ms');
  h_ylab=ylabel('Number of injections');
  grid on
  set(h_xlab,'FontSize',16,'FontName','Times');
  set(h_ylab,'FontSize',16,'FontName','Times');
  set(gca,'FontSize',16,'FontName','Times');
  saveas(gcf,'histH1H2L1dt.png')


  %%%%%%%%%%%%%%%% FREQUENCY %%%%%%%%%%%%%%%
  % df/<f> v d_eff
  figure
  eval(['semilogx(fsim.d' site '(inj' ifo 't.ind),2.*(inj' ifo 't.f-fsim.f(inj' ifo 't.ind))./(inj' ifo 't.f+fsim.f(inj' ifo 't.ind)),''.'')'])
  grid on
  xlabel('d_{eff_{inj}} / Mpc')
  ylabel('2 * (f_{0_{det}}-f_{0_{inj}}) / (f_{0_{det}}+f_{0_{inj}})')
  eval(['title(''' ifo ' Fractional difference in detected and injected central frequency versus injected effective distance'')'])
  eval(['saveas(gcf,''' ifo 'injdetfvd.png'')'])

  % df/<f> v finj 
  figure
  eval(['semilogx(fsim.f(inj' ifo 't.ind),2.*(inj' ifo 't.f-fsim.f(inj' ifo 't.ind))./(inj' ifo 't.f+fsim.f(inj' ifo 't.ind)),''.'')'])
  grid on
  xlabel('f_{0inj} / Hz')
  ylabel('2 * (f_{0_{det}}-f_{0_{inj}}) / (f_{0_{det}}+f_{0_{inj}})')
  hold on
  plot([50,50],[-2,2],'k')
  plot([2000,2000],[-2,2],'k')
  hold off
  eval(['title(''' ifo ': Fractional difference in detected and injected central frequency versus injected frequency'')'])
  eval(['saveas(gcf,''' ifo 'injdetfvf.png'')'])

  % df/<f> v log(f) v snr
  figure
  eval(['scatter(log10(fsim.f(inj' ifo 't.ind)),2.*(inj' ifo 't.f-fsim.f(inj' ifo 't.ind))./(inj' ifo 't.f+fsim.f(inj' ifo 't.ind)),30,log10(inj' ifo 't.snr),''.'')'])
  colorbar
  grid on
  xlabel('log_{10}(f_{0inj} / Hz)')
  ylabel('2 * (f_{0_{det}}-f_{0_{inj}}) / (f_{0_{det}}+f_{0_{inj}})')
  hold on
  plot([log10(50),log10(50)],[-2,2],'k')
  plot([log10(2000),log10(2000)],[-2,2],'k')
  hold off
  eval(['title(''' ifo ': Fractional difference in detected and injected central frequency versus injected frequency versus snr'')'])
  eval(['saveas(gcf,''' ifo 'injdetfvfvsnr.png'')'])

  % histogram of df/<f>
  figure
  eval(['hist(2.*(inj' ifo 't.f-fsim.f(inj' ifo 't.ind))./(inj' ifo 't.f+fsim.f(inj' ifo 't.ind)),100)'])
  grid on
  xlabel('2 * (f_{0_{det}}-f_{0_{inj}}) / (f_{0_{det}}+f_{0_{inj}})')
  eval(['title(''' ifo ': Fractional difference in detected and injected central frequency'')'])
  eval(['saveas(gcf,''' ifo 'histinjdetf.png'')'])

  %%%%%%%%%%%%% QUALITY %%%%%%%%%%%%
  % Q_det-Qinj v d_eff
  figure
  eval(['semilogx(fsim.d' site '(inj' ifo 't.ind),inj' ifo 't.Q-fsim.q(inj' ifo 't.ind),''.'')'])
  grid on
  xlabel('d_{eff_{inj}} / Mpc')
  ylabel('Q_{det}-Q_{inj}')
  eval(['title(''' ifo ': detected-injected quality factor versus injected effective distance'')'])
  eval(['saveas(gcf,''' ifo 'injdetQvd.png'')'])

  % dQ v Q_inj
  figure
  eval(['plot(fsim.q(inj' ifo 't.ind),inj' ifo 't.Q-fsim.q(inj' ifo 't.ind),''ro'')'])
  grid on
  xlabel('Q_{inj}')
  ylabel('Q_{det}-Q_{inj}')
  eval(['title(''' ifo ': detected-injected quality factor versus injected quality'')'])
  eval(['saveas(gcf,''' ifo 'injdetQvQinj.png'')'])

  % histogram of dQ
  figure
  eval(['hist(inj' ifo 't.Q-fsim.q(inj' ifo 't.ind),100)'])
  grid on
  xlabel('Q_{det}-Q_{inj}')
  eval(['title(''' ifo ': detected-injected quality factor'')'])
  eval(['saveas(gcf,''' ifo 'histinjdetQ.png'')'])

  %%%%%%%%%%%%% EFFECTIVE DISTANCE %%%%%%%%%%%%%%%
  % dd_eff/<d_eff> vs d_eff
  figure
  eval(['semilogx(fsim.d' site '(inj' ifo 't.ind),2.*(inj' ifo 't.d-fsim.d' site '(inj' ifo 't.ind))./(inj' ifo 't.d+fsim.d' site '(inj' ifo 't.ind)),''.'')'])
  grid on
  xlabel('d_{eff_{inj}} / Mpc')
  ylabel('2 * (d_{eff_{det}}-d_{eff_{inj}}) / (d_{eff_{det}}+d_{eff_{inj}})')
  eval(['title(''' ifo ': Fractional difference in detected and injected effective distance versus injected effective distance'')'])
  eval(['saveas(gcf,''' ifo 'injdetdvd.png'')'])

  % dd_eff/<d_eff> v log(f_inj)
  figure
  h_p=eval(['semilogx(fsim.f(inj' ifo 't.ind),2.*(inj' ifo 't.d-fsim.d' site '(inj' ifo 't.ind))./(inj' ifo 't.d+fsim.d' site '(inj' ifo 't.ind)),''.'')']);
  grid on
  h_xlab=xlabel('f_{inj} / Hz');
  h_ylab=ylabel('\deltad_{eff} / <d_{eff}>');
  eval(['title(''' ifo ': Fractional difference in detected and injected effective distance versus injected frequency'')'])
  grid on
  set(h_xlab,'FontSize',16,'FontName','Times');
  set(h_ylab,'FontSize',16,'FontName','Times');
  set(gca,'FontSize',16,'FontName','Times');
  set(h_p,'MarkerSize',6)
  eval(['saveas(gcf,''' ifo 'injdetdvf.png'')'])

  % hist of dd_eff/<d_eff>
  figure
  eval(['hist(2.*(inj' ifo 't.d-fsim.d' site '(inj' ifo 't.ind))./(inj' ifo 't.d+fsim.d' site '(inj' ifo 't.ind)),100)'])
  grid on
  xlabel('2 * (d_{eff_{det}}-d_{eff_{inj}}) / (d_{eff_{det}}+d_{eff_{inj}})')
  eval(['title(''' ifo ': Fractional difference in detected and injected effective distance'')'])
  eval(['saveas(gcf,''' ifo 'histinjdetd.png'')'])

  [aH1,bH1]=hist(2.*(injH1t.d-fsim.dh(injH1t.ind))./(injH1t.d+fsim.dh(injH1t.ind)),100);
  [aH2,bH2]=hist(2.*(injH2t.d-fsim.dh(injH2t.ind))./(injH2t.d+fsim.dh(injH2t.ind)),100);
  [aL1,bL1]=hist(2.*(injL1t.d-fsim.dl(injL1t.ind))./(injL1t.d+fsim.dl(injL1t.ind)),100);
  h_H1=plot(bH1,aH1,'r');
  hold on
  h_H2=plot(bH2,aH2,'b');
  h_L1=plot(bL1,aL1,'g');

  bin=-2:4/100:2;
  [aH1,bH1]=hist(2.*(injH1t.d-fsim.dh(injH1t.ind))./(injH1t.d+fsim.dh(injH1t.ind)),bin);
  [aH2,bH2]=hist(2.*(injH2t.d-fsim.dh(injH2t.ind))./(injH2t.d+fsim.dh(injH2t.ind)),bin);
  [aL1,bL1]=hist(2.*(injL1t.d-fsim.dl(injL1t.ind))./(injL1t.d+fsim.dl(injL1t.ind)),bin);

  figure
  h_H1=plot(bH1,aH1,'r');
  hold on
  h_H2=plot(bH2,aH2,'b');
  h_L1=plot(bL1,aL1,'g');
  set(h_H1,'LineWidth',2)
  set(h_H2,'LineWidth',2)
  set(h_L1,'LineWidth',2)
  h_leg=legend('H1','H2','L1');
  h_xlab=xlabel('\deltad_{eff} / <d_{eff}>');
  h_ylab=ylabel('Number of injections');
  grid on
  set(h_xlab,'FontSize',16,'FontName','Times');
  set(h_ylab,'FontSize',16,'FontName','Times');
  set(gca,'FontSize',16,'FontName','Times');
  saveas(gcf,'histH1H2L1ddd.png')

  cd ..
