function param_accuracy( veto_level, injtype, ifo, injected_list, detected_list )

%
% NULL = param_accuracy ( veto_level, injtype, ifo, injected_list, detected_list )
%
% veto_level='CAT2', 'CAT23','CAT234'
% injtype='RINGDOWN', 'EOBNR', 'PHENOM'
% The injected_list will be the sim table values.
% The detected_list will be named inj*, bg*, or pg*.
%
%
% Sarah Caudill, Oct 28th 2009

  % read in the injection list and corresponding detected list
  eval(['inj=load(''' injected_list ''');'])
  eval(['det' ifo '=load(''' detected_list ''');'])
 
  % Find number of injections in file.
  eval(['num_inj=length(det' ifo '.q);'])

 
  %============== Injected vs Detected Parameter Plots ==============%
 
  %%%%%%%%%%%%%%%%%% INJ FREQUENCY VS DET FREQ %%%%%%%%%%%%%%%%
 
  figure
  eval(['loglog(inj.f, det' ifo '.f, ''k.'')'])
  hold on
  grid on
  x_lab=xlabel('Freq_{inj} (Hz)');
  y_lab=ylabel('Freq_{det} (Hz)');
  set(x_lab,'FontSize',14);
  set(y_lab,'FontSize',14);
  set(gca,'FontSize',14);
  eval(['plot_title=title(''' veto_level ' ' injtype ': ' ifo ' Freq_{det} versus ' ifo ' Freq_{inj}'');'])
  set(plot_title,'FontSize',16,'FontWeight','b');
  eval(['saveas(gcf,''' veto_level '_' injtype '_' ifo '_detvsinj_freq.png'')'])


  %%%%%%%%%%%%%%%%%% INJ MASS VS DET MASS %%%%%%%%%%%%%%%%

  figure
  eval(['loglog(inj.m, det' ifo '.m, ''k.'')'])
  hold on
  grid on
  x_lab=xlabel('Mass_{inj} (M_{sun})');
  y_lab=ylabel('Mass_{det} (M_{sun})');
  set(x_lab,'FontSize',14);
  set(y_lab,'FontSize',14);
  set(gca,'FontSize',14);
  eval(['plot_title=title(''' veto_level ' ' injtype ': ' ifo ' Mass_{det} versus ' ifo ' Mass_{inj}'');'])
  set(plot_title,'FontSize',16,'FontWeight','b');
  eval(['saveas(gcf,''' veto_level '_' injtype '_' ifo '_detvsinj_mass.png'')'])
 
 
  %%%%%%%%%%%%%%%%%% TIME HISTOGRAM %%%%%%%%%%%%%%%%
 
  if(strcmp(ifo,'H1') || strcmp(ifo,'H2'))
    figure
    eval(['hist(inj.th - det' ifo '.t, 100)'])
    hold on
    grid on
    x_lab=xlabel('Time_{inj} - Time_{det} (s)');
    y_lab=ylabel('Number');
    set(x_lab,'FontSize',14);
    set( y_lab,'FontSize',14);
    set(gca,'FontSize',14);
    eval(['plot_title=title(''' veto_level ' ' injtype ': Histogram of ' ifo ' Time_{inj} - Time_{det}'');'])
    set(plot_title,'FontSize',16,'FontWeight','b');
    eval(['saveas(gcf,''' veto_level '_' injtype '_' ifo '_injminusdet_timehist.png'')'])
  else
    figure
    eval(['hist(inj.tl - det' ifo '.t, 100)'])
    hold on
    grid on
    x_lab=xlabel('Time_{inj} - Time_{det} (s)');
    y_lab=ylabel('Number');
    set(x_lab,'FontSize',14);
    set(y_lab,'FontSize',14);
    set(gca,'FontSize',14);
    eval(['plot_title=title(''' veto_level ' ' injtype ': Histogram of ' ifo ' Time_{inj} - Time_{det}'');'])
    set(plot_title,'FontSize',16,'FontWeight','b');
    eval(['saveas(gcf,''' veto_level '_' injtype '_' ifo '_injminusdet_timehist.png'')'])
  end
 
 
  %%%%%%%%%%%%%%%%%% TIME VS FREQUENCY %%%%%%%%%%%%%%%%
 
  if(strcmp(ifo,'H1') || strcmp(ifo,'H2'))
    figure
    eval(['semilogx(inj.f, inj.th - det' ifo '.t, ''k.'')'])
    hold on
    grid on
    x_lab=xlabel('Freq_{inj} (Hz)');
    y_lab=ylabel('Time_{inj} - Time_{det} (s)');
    set(x_lab,'FontSize',14);
    set(y_lab,'FontSize',14);
    set(gca,'FontSize',14);
    eval(['plot_title=title(''' veto_level ' ' injtype ': ' ifo ' Time_{inj} - Time_{det} versus ' ifo ' Freq_{inj}'');'])
    set(plot_title,'FontSize',16,'FontWeight','b');
    eval(['saveas(gcf,''' veto_level '_' injtype '_' ifo '_detvsinj_timefreq.png'')'])
  else
    figure
    eval(['semilogx(inj.f, inj.tl - det' ifo '.t, ''k.'')'])
    hold on
    grid on
    x_lab=xlabel('Freq_{inj} (Hz)');
    y_lab=ylabel('Time_{inj} - Time_{det} (s)');
    set(x_lab,'FontSize',14);
    set(y_lab,'FontSize',14);
    set(gca,'FontSize',14);
    eval(['plot_title=title(''' veto_level ' ' injtype ': ' ifo ' Time_{inj} - Time_{det} versus ' ifo ' Freq_{inj}'');'])
    set(plot_title,'FontSize',16,'FontWeight','b');
    eval(['saveas(gcf,''' veto_level '_' injtype '_' ifo '_detvsinj_timefreq.png'')'])
  end
 
 
  %%%%%%%%%%%%%%%%%% INJ EFF DISTANCE VS DET EFF DISTANCE %%%%%%%%%%%%%%%%
 
  if(strcmp(ifo,'H1') || strcmp(ifo,'H2'))
    figure
    eval(['loglog(inj.dh, det' ifo '.d, ''k.'')'])
    hold on
    grid on
    xpts=logspace(-2,3,500);
    loglog(xpts, xpts,'r-');
    x_lab=xlabel('Eff Distance_{inj} (Mpc)');
    y_lab=ylabel('Eff Distance_{det} (Mpc)');
    set(x_lab,'FontSize',14);
    set(y_lab,'FontSize',14);
    set(gca,'FontSize',14);
    eval(['plot_title=title(''' veto_level ' ' injtype ': ' ifo ' Eff Distance_{det} versus ' ifo ' Eff Distance_{inj}'');'])
    set(plot_title,'FontSize',16,'FontWeight','b');
    eval(['saveas(gcf,''' veto_level '_' injtype '_' ifo '_detvsinj_distance.png'')'])
  else
    figure
    eval(['loglog(inj.dl, det' ifo '.d, ''k.'')'])
    hold on
    grid on
    xpts=logspace(-2,3,500);
    loglog(xpts, xpts,'r-');
    x_lab=xlabel('Eff Distance_{inj} (Mpc)');
    y_lab=ylabel('Eff Distance_{det} (Mpc)');
    set(x_lab,'FontSize',14);
    set(y_lab,'FontSize',14);
    set(gca,'FontSize',14);
    eval(['plot_title=title(''' veto_level ' ' injtype ': ' ifo ' Eff Distance_{det} versus ' ifo ' Eff Distance_{inj}'');'])
    set(plot_title,'FontSize',16,'FontWeight','b');
    eval(['saveas(gcf,''' veto_level '_' injtype '_' ifo '_detvsinj_distance.png'')'])
  end
 
 
  %%%%%%%%%%%%%%%%%% FRACTIONAL EFF DISTANCE %%%%%%%%%%%%%%%%
 

  if(strcmp(ifo,'H1') || strcmp(ifo,'H2'))
    figure
    eval(['semilogx(inj.dh, (inj.dh-det' ifo '.d)./((inj.dh+det' ifo '.d)./2), ''k.'')'])
    hold on
    grid on
    x_lab=xlabel('Eff Distance_{inj} (Mpc)');
    y_lab=ylabel('(Eff Distance_{inj} - Eff Distance_{det}) / <Eff Distance>');
    set(x_lab,'FontSize',14);
    set(y_lab,'FontSize',14);
    set(gca,'FontSize',14);
    eval(['plot_title=title(''' veto_level ' ' injtype ': Fractional Effective Distance Accuracy in ' ifo ''');'])
    set(plot_title,'FontSize',16,'FontWeight','b');
    eval(['saveas(gcf,''' veto_level '_' injtype '_' ifo '_accuracyvsinj_distance.png'')'])
  else
    figure
    eval(['semilogx(inj.dl, (inj.dl-det' ifo '.d)./((inj.dl+det' ifo '.d)./2), ''k.'')'])
    hold on
    grid on
    x_lab=xlabel('Eff Distance_{inj} (Mpc)');
    y_lab=ylabel('(Eff Distance_{inj} - Eff Distance_{det}) / <Eff Distance>');
    set(x_lab,'FontSize',14);
    set(y_lab,'FontSize',14);
    set(gca,'FontSize',14);
    eval(['plot_title=title(''' veto_level ' ' injtype ': Fractional Effective Distance Accuracy in ' ifo ''');'])
    set(plot_title,'FontSize',16,'FontWeight','b');
    eval(['saveas(gcf,''' veto_level '_' injtype '_' ifo '_accuracyvsinj_distance.png'')'])
  end
 
 
  %%%%%%%%%%%%%%%%%% FRACTIONAL FREQUENCY %%%%%%%%%%%%%%%%
 
  figure
  eval(['semilogx(inj.f, (inj.f-det' ifo '.f)./((inj.f+det' ifo '.f)./2), ''k.'')'])
  hold on
  grid on
  x_lab=xlabel('Freq_{inj} (Hz)');
  y_lab=ylabel('(Freq_{inj} - Freq_{det}) / <Freq>');
  set(x_lab,'FontSize',14);
  set(y_lab,'FontSize',14);
  set(gca,'FontSize',14);
  eval(['plot_title=title(''' veto_level ' ' injtype ': Fractional Frequency Accuracy in ' ifo ''');'])
  set(plot_title,'FontSize',16,'FontWeight','b');
  eval(['saveas(gcf,''' veto_level '_' injtype '_' ifo '_accuracyvsinj_frequency.png'')'])
