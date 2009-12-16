function doub_hist( veto_level, coinctype, ifotime, ifo1, ifo2, ifo1bgfile, ifo2bgfile, nonbgtype, ifo1nonbgfile, ifo2nonbgfile )

%
% veto_level = 'CAT2','CAT23', 'CAT234'
%
% coinctype = 'doub', 'trip'
%
% ifotime = 'doub', 'trip'
%
% ifo1, ifo2 = 'H1', 'H2', 'L1'
%
% ifo1bgfile, etc. are the background *doub.mat files for each ifo.
%
% nonbgtype = 'pg','fg'
%
% ifo1nonbgfile, etc. are the playground or foreground *doub.mat files for each ifo.
%
% Sarah Caudill, Nov 16th 2009

  % read in the background and playground or foreground files
  eval(['ifo1bg=load(''' ifo1bgfile ''');'])
  eval(['ifo2bg=load(''' ifo2bgfile ''');'])
  eval(['ifo1nonbg=load(''' ifo1nonbgfile ''');'])
  eval(['ifo2nonbg=load(''' ifo2nonbgfile ''');'])

  % call function which counts the number of events per timeslide and zero lag and makes a histogram
  eval(['plothist(veto_level, coinctype, ifotime, ifo1, ifo2, ''' ifo1bgfile ''', nonbgtype, ''' ifo1nonbgfile ''');'])

  % calculate the ranking for background based on triple det stat
  bgdetstat = min(2.*ifo1bg.snr + 2.2, (min(2.*ifo2bg.snr + 2.2, ifo1bg.snr + ifo2bg.snr)));

  % calculate the ranking for pg/zerolag based on triple det stat
  nonbgdetstat = min(2.*ifo1nonbg.snr + 2.2, (min(2.*ifo2nonbg.snr + 2.2, ifo1nonbg.snr + ifo2nonbg.snr)));

  % call another function that plots cumulative hist
  plotcumhist(veto_level, coinctype, ifotime, ifo1, ifo2, bgdetstat, nonbgtype, nonbgdetstat);

  % call another function that plots DT vs F plots for timeslides
  eval(['background_param_plots(veto_level, coinctype, ifotime, ifo1, ifo2,'''  ifo1bgfile ''',''' ifo2bgfile ''');'])
