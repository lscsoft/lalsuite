function background_param_plots( veto_level,coinctype, ifotime, ifo1, ifo2, bgfile1, bgfile2 )

%
% veto_level = 'CAT2','CAT23','CAT234'
%
% coinctype = 'trip', 'doub'
%
% ifotime = 'trip', 'doub'
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
  eval(['plot_title=title([''dt (' ifo1 '' ifo2 ') v. f for ' coinctype ' coincidence background events in ' ifotime ' time'']);'])
  set(plot_title,'FontSize',14,'FontWeight','b');
  grid on;
  if(strcmp(coinctype,'trip') && strcmp(ifotime,'trip'))
    eval(['saveas(gcf,''' veto_level '_' ifo1 '' ifo2 '_H1H2L1_H1H2L1_bg_dtvsf.png'')'])
  elseif(strcmp(coinctype,'doub') && strcmp(ifotime,'trip'))
    eval(['saveas(gcf,''' veto_level '_' ifo1 '' ifo2 '_' ifo1 '' ifo2 '_H1H2L1_bg_dtvsf.png'')'])
  else
    eval(['saveas(gcf,''' veto_level '_' ifo1 '' ifo2 '_' ifo1 '' ifo2 '_' ifo1 '' ifo2 '_bg_dtvsf.png'')'])
  end
  close;
