function plotcumhist( veto_level, coinctype, ifotime, ifo1, ifo2, bgdetstat, nonbgtype, nonbgdetstat )

%
% null = plotcumhist( veto_level, coinctype, ifotime, ifo1, ifo2, bgdetstat, nonbgdetstat )
%
% veto_level = 'CAT2','CAT23','CAT234'
%
% coinctype = 'trip', 'doub'
%
% ifo1, ifo2 can both be 'trip' if plothist is passed its options
% through trip_hist or one of 'H1', 'H2', or 'L1' if plothist
% is passed its options through doub_hist.
%
% bgdetstat is the detection statistic calculated in trip_hist or doub_hist
% for background triggers.
%
% nonbgtype = 'pg','fg'
%
% nonbgdetstat is the detection statistic calculated in trip_hist or doub_hist
% for zerolag or playground triggers.
%
% Sarah Caudill, Nov 16th 2009

  %%%%%%%%%%%%% CUMULATIVE HISTOGRAMS OF ZERO-LAG AND BACKGROUND COINC %%%%%%%%%%%%%

  %plotting range
  limit=max(max(bgdetstat),max(nonbgdetstat));
  % binning the det stat
  rho=0:limit/50:limit;

  for i=1:length(rho)
    % count the number of coincs below a given detstat value
    Nbg(i)=length(bgdetstat(bgdetstat>rho(i)));
  end
  %scaling the background to one "experiment"
  Nbg=Nbg./100;

  for i=1:length(rho)
    Nnonbg(i)=length(nonbgdetstat(nonbgdetstat>rho(i)));
  end
  %scaling playground histogram to full data
  if(strcmp(nonbgtype,'pg'))
    Nnonbg(i)=length(nonbgdetstat(nonbgdetstat>rho(i))).*6370./600;
  end

  for i=1:length(Nnonbg)
    if Nnonbg(i)==0;
       % want to plot on a log scale, and so have to set values that are =0 to a small non-xero number
       Nnonbg(i)=1e-5;
    end
  end

  for i=1:length(Nbg)
    if Nbg(i)==0;
      Nbg(i)=1e-5;
    end
  end


  % error bars
  errp=Nbg+sqrt(Nbg);
  errm=Nbg-sqrt(Nbg);
  % scaling playground errorbars to full data error bars
  if(strcmp(nonbgtype,'pg'))
    errp=Nbg+sqrt(Nbg).*6370./600;
    errm=Nbg-sqrt(Nbg).*6370./600;
  end

  for i=1:length(errm)
    if errm(i)<=0;
      errm(i)=1e-5;
    end
  end

  figure %make the cumulative histogram
  h_bg=semilogy(rho,Nbg,'kx');
  grid on
  hold on
  h_nonbg=semilogy(rho,Nnonbg,'ro');
  for i=1:length(rho)
    h_er=semilogy([rho(i) rho(i)],[errm(i) errp(i)],'g.-');
    set(h_er,'MarkerSize',15,'LineWidth',2);
  end
  set(h_bg,'MarkerSize',15);
  set(h_nonbg,'MarkerSize',10);
  h_xlab=xlabel('\rho_{det-stat}');
  h_ylab=ylabel('Number of Coincidences');
  % title('Number of triggers above detection statistic, H1L1 background and intime doubles in triple time')
  if(strcmp(nonbgtype,'pg'))
    h_leg=legend('background/100','playground zero-lag','1\sigma');
  else
    h_leg=legend('background/100','zero-lag','1\sigma');
  end
  set(h_xlab,'FontSize',14);
  set(h_ylab,'FontSize',14);
  set(h_leg,'FontSize',14);
  set(gca,'FontSize',14);
  axis([0 1 1e-3 1e3]);
  axis autox
 
  if(strcmp(coinctype,'trip') && strcmp(ifotime,'trip'))
    eval(['saveas(gcf,''' veto_level '_H1H2L1_H1H2L1_bg' nonbgtype '_Ndetstat.png'')'])
  elseif(strcmp(coinctype,'doub') && strcmp(ifotime,'trip'))
    eval(['saveas(gcf,''' veto_level '_' ifo1 '' ifo2 '_H1H2L1_bg' nonbgtype '_Ndetstat.png'')'])
  else
    eval(['saveas(gcf,''' veto_level '_' ifo1 '' ifo2 '_' ifo1 '' ifo2 '_bg' nonbgtype '_Ndetstat.png'')'])
  end
  close;
