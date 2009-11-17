function scatter_plots( veto_level, coinctype, ifo1, ifo2, bgfile1, bgfile2, injtype, nonbgtype, nonbgfile1, nonbgfile2 )

%
% NULL = scatter_plots ( veto_level, coinctype, ifo1, ifo2, bgfile1, bgfile2, injtype,  nonbgtype, nonbgfile1, nonbgfile2 )
%
% coinctype = 'trip', 'doub'
%
% veto_level = 'NOVETO','CAT2','CAT23'
%
% ifo1='H1','H2','L1'
%
% bgfile1 and bgfile2 should always be background.mat files. They will be plotted
% first. bgfile1 will be plotted on the x-axis and bgfile2 will be plotted on
% y-axis.
%
% injtype = 'RINGDOWN', 'EOBNR', 'PHENOM', 'NONE'
%
% nonbgtype = 'inj','pg','fg'
%
% nonbgfile1 will be plotted on x-axis of plots. nonbgfile2 will be plotted on y-axis.
% These are non-background files. They will be plotted second so they will be injection
% files, playground files, or foreground files.
%
% Sarah Caudill, Nov 16th 2009

  % read in the background and injection or playground files
  eval(['bgx=load(''' bgfile1 ''');'])
  eval(['bgy=load(''' bgfile2 ''');'])
  eval(['nonbgx=load(''' nonbgfile1 ''');'])
  eval(['nonbgy=load(''' nonbgfile2 ''');'])


  %======== Timeslides and Injection Plots ========%

    %%%%%%%%%%%%%%%%%% SNR vs SNR SCATTERPLOTS %%%%%%%%%%%%%%%%

  if(strcmp(nonbgtype,'inj'))
    figure
    loglog(bgx.snr, bgy.snr, 'kx')
    hold on
    loglog(nonbgx.snr, nonbgy.snr, 'rx')
    grid on
    eval(['x_lab=xlabel(''\rho_{' ifo1 '}'');'])
    eval(['y_lab=ylabel(''\rho_{' ifo2 '}'');'])
    set(x_lab,'FontSize',14);
    set(y_lab,'FontSize',14);
    set(gca,'FontSize',14);
    eval(['plot_title=title(''' veto_level ' ' injtype ': ' ifo2 ' SNR versus ' ifo1 ' SNR'');'])
    set(plot_title,'FontSize',14,'FontWeight','b');
    eval(['leg=legend(''Background ' coinctype ''',''Injections ' coinctype ''',0);'])
    set(leg,'FontSize',14);
    eval(['saveas(gcf,''' veto_level '_' injtype '_' ifo1 '' ifo2 '_' coinctype '_' nonbgtype '_snrscatterplot.png'')'])
    close;

  %======== Timeslides and Playground or Foreground Plots ========%

    %%%%%%%%%%%%%%%%%% SNR vs SNR SCATTERPLOTS %%%%%%%%%%%%%%%%

  else
    figure
    loglog(bgx.snr, bgy.snr, 'kx')
    hold on
    loglog(nonbgx.snr, nonbgy.snr, 'bx')
    grid on
    eval(['x_lab=xlabel(''\rho_{' ifo1 '}'');'])
    eval(['y_lab=ylabel(''\rho_{' ifo2 '}'');'])
    set(x_lab,'FontSize',14);
    set(y_lab,'FontSize',14);
    set(gca,'FontSize',14);
    eval(['plot_title=title(''' veto_level ': ' ifo2 ' SNR versus ' ifo1 ' SNR'');'])
    set(plot_title,'FontSize',14,'FontWeight','b');
    if(strcmp(nonbgtype,'pg'))
      eval(['leg=legend(''Background ' coinctype ''',''Playground ' coinctype ''',0);'])
    elseif(strcmp(nonbgtype,'fg'))
      eval(['leg=legend(''Background ' coinctype ''',''Zero-lag ' coinctype ''',0);'])
    else
      eval(['leg=legend(''Background ' coinctype ''',''' nonbgtype ' ' coinctype ''',0);'])
    end
    set(leg,'FontSize',14);
    eval(['saveas(gcf,''' veto_level '_' injtype '_' ifo1 '' ifo2 '_' coinctype '_' nonbgtype '_snrscatterplot.png'')'])
    close;
  end
