function background_param_plots( veto_level,coinctype, ifo1, ifo2, bgfile1, bgfile2 )

%
% veto_level = 'NOVETO','CAT2','CAT23'
%
% coinctype = 'trip', 'doub'
%
% ifo1, ifo2 = 'H1', 'H2', 'L1'
%
% bgfile1, bgfile2 are the background *.mat files.
%
%% Sarah Caudill, Nov 16th 2009


  % read in the background files
  eval(['bg1=load(''' bgfile1 ''');'])
  eval(['bg2=load(''' bgfile2 ''');'])

  %%%%%%%%%%%%% DT VS F PLOTS %%%%%%%%%%%%%

  % plot dt v. f for all events
  dt = mod( (bg1.t-bg2.t), 1 ); % reslide data
  dt = (dt - (dt > .5)).*1.e3; % adjust for t2>t1 and convert to s -> ms
  figure
  semilogx(bg1.f,dt,'b+');
  x_lab=xlabel('f (Hz)');
  eval(['y_lab=ylabel(''dt between ' ifo1 ' & ' ifo2 ' (ms)'');'])
  set(x_lab,'FontSize',14);
  set(y_lab,'FontSize',14);
  set(gca,'FontSize',14);
  eval(['plot_title=title([''dt (' ifo1 '' ifo2 ') v. f for ' coinctype ' time background events'']);'])
  set(plot_title,'FontSize',14,'FontWeight','b');
  grid on;
  eval(['saveas(gcf,''' veto_level '_' ifo1 '' ifo2 '_' coinctype 'bgdtvsf.png'')'])
  close;
