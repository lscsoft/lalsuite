
function ifoifo( found_list, ifo1, ifo2 )

%
% NULL = ifoifo ( found_list, ifo1, ifo2 )
%
% Make parameter estimation plots comparing values recovered in two
% different ifo's.  found_list is a column vector af filenames output
% from lalapps_coincringread, ifo1 and ifo2 are two character strings
% for the ifos you want to compare.  Certain three ifo plots are
% always generated.
%
% Plots are saved in paramacc, which is generated if it does not exist.
%
% EXAMPLE:
% 
%    found_list = [ 'injH1H2L1coincs_1001.xml'; 'injH1H2L1coincs_1002' ];
%    ifoifo(found_list, 'H1', 'H2');
%
% Lisa Goggin, Mar 15th 2006
% modified: P. T. Baker, 7 July 2009

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
  
  %=========================== Injected Quantities ===========================%

  % Read found sim_ringdown tables
  N_files = length(found_list(:,1));

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

  fsim.ind=1:length(fsim.th);

  % transpose fsim arrays                       
  fsim.th=transpose(fsim.th);
  fsim.tl=transpose(fsim.tl);
  fsim.f=transpose(fsim.f);
  fsim.Q=transpose(fsim.q);
  fsim.dh=transpose(fsim.dh);
  fsim.dl=transpose(fsim.dl);
  fsim.ind=transpose(fsim.ind);

  % get sites from input ifos
  if strcmp(ifo1,'H1') || strcmp(ifo1,'H2')
    site1 = 'h';
  elseif strcmp(ifo1,'L1')
    site1 = 'l';
  end

  if strcmp(ifo2,'H1') || strcmp(ifo2,'H2')
    site2 = 'h';
  elseif strcmp(ifo2,'L1')
    site2 = 'l';
  end

  % ds^2
  eval(['injdetds2=get_ds_sq(inj' ifo1 't.f,inj' ifo1 't.Q,inj' ifo2 't.f,inj' ifo2 't.Q);'])

  eval(['inj' ifo1 't.tau=inj' ifo1 't.Q./pi./inj' ifo1 't.f;'])
  eval(['inj' ifo2 't.tau=inj' ifo2 't.Q./pi./inj' ifo2 't.f;'])

  if ~exist('paramacc','dir')
    mkdir paramacc
  end
  cd paramacc

  %%%%%%%%%%%%%%%%%%%%%%% difference scatter plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%  TIME %%%%%%%%%%%%%%%%%%%%%%
  % dt vs distance
  figure
  eval(['semilogx(inj' ifo1 't.d,inj' ifo1 't.t-inj' ifo2 't.t,''.'')'])
  grid on
  eval(['xlabel(''d_{eff_{' ifo1 '}} / Mpc'')'])
  eval(['ylabel(''t_{' ifo1 '}-t_{' ifo2 '} / s'')'])
  eval(['title(''' ifo1 '-' ifo2 ' time versus ' ifo1 ' effective distance'')'])
  eval(['saveas(gcf,''' ifo1 '' ifo2 'ttvd.png'')'])

  % dt vs tau
  figure
  eval(['semilogx(inj' ifo1 't.tau,inj' ifo1 't.t-inj' ifo2 't.t,''.'')'])
  grid on
  eval(['xlabel(''\tau_{' ifo1 '} / s'')'])
  eval(['ylabel(''(t_{' ifo1 '}-t_{' ifo2 '}) / s'')'])
  eval(['title(''' ifo1 '-' ifo2 ' time versus ' ifo1 ' \tau'')'])
  eval(['saveas(gcf,''' ifo1 '' ifo2 'ttvtau.png'')'])

  % dt (in ms) vs f
  figure
  eval(['semilogx(inj' ifo1 't.f,1000*(inj' ifo1 't.t-inj' ifo2 't.t),''.'')'])
  grid on
  h_xlab=eval(['xlabel(''f_{' ifo1 '} / Hz'')']);
  %h_ylab=eval(['ylabel(''(t_{' ifo1 '}-t_{' ifo2 '}) / ms'')']);
  h_ylab=ylabel('\delta t / ms');
  %set(h_xlab,'FontSize',16,'FontName','Times');
  %set(h_ylab,'FontSize',16,'FontName','Times');
  %set(gca,'FontSize',16,'FontName','Times');
  eval(['title(''' ifo1 '-' ifo2 ' time versus ' ifo1 ' frequency'')'])
  eval(['saveas(gcf,''' ifo1 '' ifo2 'ttvf.png'')'])

  figure
  semilogx(injH1t.f,1000*(injH1t.t-injL1t.t),'r.')
  hold on
  semilogx(injH2t.f,1000*(injH2t.t-injL1t.t),'b.')
  grid on
  h_xlab=xlabel('f_{H1} / Hz');
  %h_ylab=eval(['ylabel(''(t_{' ifo1 '}-t_{' ifo2 '}) / ms'')']);
  h_ylab=ylabel('\delta t / ms');
  h_leg=legend('H1-L1','H2-L1');
  %set(h_xlab,'FontSize',16,'FontName','Times');
  %set(h_ylab,'FontSize',16,'FontName','Times');
  %set(gca,'FontSize',16,'FontName','Times');
  eval(['title(''H1-L1 and H2-L1  time versus ' ifo1 ' frequency'')'])
  saveas(gcf,'H1L1H2L1ttvf.png')


  %%%%%%%%%%%%%%%%%% FREQUENCY %%%%%%%%%%%%%%%%
  % df/<f> vs d
  figure
  eval(['semilogx(inj' ifo1 't.d,2.*(inj' ifo1 't.f-inj' ifo2 't.f)./(inj' ifo1 't.f+inj' ifo2 't.f),''.'')'])
  grid on
  eval(['xlabel(''d_{eff_{' ifo1 '}} / Mpc'')'])
  eval(['ylabel(''2 * (f_{0_{' ifo1 '}}-f_{0_{' ifo2 '}}) / (f_{0_{' ifo1 '}}+f_{0_{' ifo2 '}})'')'])
  eval(['title(''Fractional difference in ' ifo1 '' ifo2 ' central frequency versus ' ifo1 ' effective distance'')'])
  eval(['saveas(gcf,''' ifo1 '' ifo2 'fvd.png'')'])

  % df/<f> vs f
  figure
  eval(['semilogx(inj' ifo1 't.f,2.*(inj' ifo1 't.f-inj' ifo2 't.f)./(inj' ifo1 't.f+inj' ifo2 't.f),''.'')'])
  grid on
  eval(['xlabel(''f_{0_{' ifo1 '}} / Hz'')'])
  eval(['ylabel(''2 * (f_{0_{' ifo1 '}}-f_{0_{' ifo2 '}}) / (f_{0_{' ifo1 '}}+f_{0_{' ifo2 '}})'')'])
  eval(['title(''Fractional difference in ' ifo1 '' ifo2 ' central frequency versus ' ifo1 ' central frequency'')'])
  eval(['saveas(gcf,''' ifo1 '' ifo2 'ffvf.png'')'])

  % df/<f> vs d vs f
  figure
  eval(['scatter(inj' ifo1 't.f,inj' ifo1 't.d,70,2.*(inj' ifo1 't.f-inj' ifo2 't.f)./(inj' ifo1 't.f+inj' ifo2 't.f),''.'')'])
  set(gca,'xscale','log')
  set(gca,'yscale','log')
  colorbar
  h_xlab=eval(['xlabel(''f_{0} / Hz'')']);
  h_ylab=eval(['ylabel(''d_{eff} / Mpc'')']);
  grid on
  %set(h_xlab,'FontSize',16,'FontName','Times'); 
  %set(h_ylab,'FontSize',16,'FontName','Times');
  %set(gca,'FontSize',16,'FontName','Times');
  eval(['title(''Fractional difference in ' ifo1 '' ifo2 ' central frequency versus ' ifo1 ' d_{eff} and central frequency'')'])
  eval(['saveas(gcf,''' ifo1 '' ifo2 'fvdvff.png'')'])


  %%%%%%%%%%%%%%%%%%% QUALITY %%%%%%%%%%%%%%%%%%%%
  % dQ vs d
  figure
  eval(['semilogx(inj' ifo1 't.d,inj' ifo1 't.Q-inj' ifo2 't.Q,''.'')'])
  grid on
  eval(['xlabel(''d_{eff_{' ifo1 '}} / Mpc'')'])
  eval(['ylabel(''Q_{' ifo1 '}-Q_{' ifo2 '}'')'])
  eval(['title(''' ifo1 '-' ifo2 ' quality factor versus ' ifo1 ' effective distance'')'])
  eval(['saveas(gcf,''' ifo1 '' ifo2 'QQvd.png'')'])


  %%%%%%%%%%%%%%%%%% EFFECTIVE DISTANCE %%%%%%%%%%%%%%%
  % dd/<d> vs d
  figure
  eval(['semilogx(inj' ifo1 't.d,2.*(inj' ifo1 't.d-inj' ifo2 't.d)./(inj' ifo1 't.d+inj' ifo2 't.d),''.'')'])
  grid on
  eval(['xlabel(''d_{eff_{' ifo1 '}} / Mpc'')'])
  eval(['ylabel(''2 * (d_{eff_{' ifo1 '}}-d_{eff_{' ifo2 '}}) / (d_{eff_{' ifo1 '}}+d_{eff_{' ifo2 '}})'')'])
  eval(['title(''Fractional difference in ' ifo1 '' ifo2 ' effective distance versus ' ifo1 ' effective distance'')'])
  eval(['saveas(gcf,''' ifo1 '' ifo2 'ddvd.png'')'])

  % dd/<d> vs f
  figure
  eval(['semilogx(inj' ifo1 't.f,2.*(inj' ifo1 't.d-inj' ifo2 't.d)./(inj' ifo1 't.d+inj' ifo2 't.d),''.'')'])
  grid on
  eval(['xlabel(''f_{' ifo1 '} / Hz'')'])
  eval(['ylabel(''2 * (d_{eff_{' ifo1 '}}-d_{eff_{' ifo2 '}}) / (d_{eff_{' ifo1 '}}+d_{eff_{' ifo2 '}})'')'])
  eval(['title(''Fractional difference in ' ifo1 '' ifo2 ' effective distance versus ' ifo1 ' frequency'')'])
  eval(['saveas(gcf,''' ifo1 '' ifo2 'ddvf.png'')'])

  % log10(|dd|) vs d vs f
  figure
  eval(['scatter(inj' ifo1 't.f,inj' ifo1 't.d,70,log10(abs(inj' ifo1 't.d-inj' ifo2 't.d)),''.'')'])
  set(gca,'xscale','log')
  set(gca,'yscale','log')
  colorbar
  h_xlab=eval(['xlabel(''f_{0} / Hz'')']);
  h_ylab=eval(['ylabel(''d_{eff} / Mpc'')']);
  grid on
  %set(h_xlab,'FontSize',20,'FontName','Times');
  %set(h_ylab,'FontSize',20,'FontName','Times');
  %set(gca,'FontSize',16,'FontName','Times');
  eval(['title(''log_{10} difference between ' ifo1 ' and ' ifo2 ' d_{eff} versus ' ifo1 ' d_{eff} and central frequency'')'])
  eval(['saveas(gcf,''' ifo1 '' ifo2 'ddvdvf.png'')'])

  cd ..
