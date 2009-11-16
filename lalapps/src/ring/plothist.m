function plothist( veto_level, coinctype, ifo1, ifo2, bgfile, nonbgtype, nonbgfile )

%
% veto_level = 'NOVETO','CAT2','CAT23'
%
% coinctype = 'trip', 'doub'
%
% ifo1, ifo2 can both be 'trip' if plothist is passed its options
% through trip_hist or one of 'H1', 'H2', or 'L1' if plothist
% is passed its options through doub_hist.
%
% bgfile is the background *.mat file.
%
% nonbgtype = 'pg','fg'
%
% nonbgfile is the playground or foreground *.mat file.
%
% Sarah Caudill, Nov 16th 2009

  %%%%%%%%%%%%% TIMESLIDE AND ZEROLAG NUMBER OF EVENTS HISTOGRAMS  %%%%%%%%%%%%%

  % read in the background and playground or foreground files
  eval(['bg=load(''' bgfile ''');'])
  eval(['nonbg=load(''' nonbgfile ''');'])
 
  % Sort the timeslides
  for k=1:length(bg.id)
    tmp=bg.id(k);
    alsotmp=tmp{1};
    bg.slidenum(k)=str2num(alsotmp(33:36));
  end

  for k=1:length(bg.slidenum)
    if bg.slidenum(k) > 5000
      event(k)=(bg.slidenum(k) - 5000)*(-1);
    else
      event(k)=bg.slidenum(k);
    end
  end

  % plot number events per slide
  figure
  hist(event,101)
  hold on
  h=findobj(gca,'Type','patch');
  % set the time-slides to red
  set(h,'FaceColor','r')
  [A,B]=hist(event,101);

  % plot number events in zerolag
  if(strcmp(nonbgtype,'pg'))
    x=0.0;
    y=length(nonbg.snr).*6370./600;
    bar(x,length(nonbg.snr).*6370./600);
    err=sqrt(length(nonbg.snr)).*6370./600;
    s=errorbar(x,y,err,'g');
    set(s,'LineStyle','none')
    set(s,'LineWidth',1)
  else
    bar(0,length(nonbg.snr));
  end

  % plot line for average number in each timeslide
  ave=sum(A)/100;
  plot([-60,60],[ave,ave],'k--')

  h_xlab=xlabel('SlideNumber');
  h_ylab=ylabel('Number of Coincidences');
  if(strcmp(ifo1,'trip'))
    eval(['h_t=title( ''Number of H1H2L1 background events in H1H2L1 time'' );'])
  else
    eval(['h_t=title( ''Number of ' ifo1 '' ifo2 '  background events in ' ifo1 '' ifo2 ' ' coinctype ' time'' );'])
  end
  set(h_xlab,'FontSize',14);
  set(h_ylab,'FontSize',14);
  set(h_t,'FontSize',14,'FontWeight','b');
  grid on
  axis([-60 60 0 320])
  set(gca,'FontSize',14);
  axis autoy
  if(strcmp(ifo1,'trip'))
    eval(['saveas(gcf,''' veto_level '_H1H2L1_' coinctype '_bg' nonbgtype 'hist.png'')'])
    close;
  else
    eval(['saveas(gcf,''' veto_level '_' ifo1 '' ifo2 '_' coinctype '_bg' nonbgtype 'hist.png'')'])
    close;
  end
  ave_num_bg = ave;

  %%%%%%%%%%%%% TIMESLIDE NUMBER OF EVENTS AS FUNCTION OF FREQ %%%%%%%%%%%%%

  figure
  hist(log10(bg.f),100)
  h_xlab=xlabel('log_{10}(f)');
  h_ylab=ylabel('Number of Coincidences');
  if(strcmp(ifo1,'trip'))
    eval(['h_t=title(''Number of H1H2L1 background events in H1H2L1 time as a function of frequency'');'])
  else
    eval(['h_t=title(''Number of ' ifo1 '' ifo2 '   background events in ' ifo1 '' ifo2 ' ' coinctype ' time as a function of frequency'');'])
  end
  set(gca,'FontSize',14);
  grid on
  set(h_xlab,'FontSize',14);
  set(h_ylab,'FontSize',14);
  set(h_t,'FontSize',14,'FontWeight','b');
  if(strcmp(ifo1, 'trip'))
    eval(['saveas(gcf,''' veto_level '_H1H2L1_' coinctype '_bg' nonbgtype 'fhist.png'')'])
    close;
  else
    eval(['saveas(gcf,''' veto_level '_' ifo1 '' ifo2 '_' coinctype '_bg' nonbgtype 'fhist.png'')'])
    close;
  end
