function trip_hist( veto_level, H1bgfile, H2bgfile, L1bgfile, nonbgtype, H1nonbgfile, H2nonbgfile, L1nonbgfile )

%
% veto_level = 'CAT2','CAT23','CAT234'
%
% H1bgfile, etc. are the background *trip.mat files for each ifo.
%
% nonbgtype = 'pg','fg'
%
% H1nonbgfile, etc. are the playground or foreground *trip.mat files for each ifo.
%
% Sarah Caudill, Nov 16th 2009

  coinctype = 'trip';
  ifotime = 'trip';
  ifo1 = 'trip';
  ifo2 = 'trip';

  % read in the background and playground or foreground files
  eval(['H1bg=load(''' H1bgfile ''');'])
  eval(['H2bg=load(''' H2bgfile ''');'])
  eval(['L1bg=load(''' L1bgfile ''');'])
  eval(['H1nonbg=load(''' H1nonbgfile ''');'])
  eval(['H2nonbg=load(''' H2nonbgfile ''');'])
  eval(['L1nonbg=load(''' L1nonbgfile ''');'])

  % call function which counts the number of events per timeslide and zero lag and makes a histogram
  eval(['plothist(veto_level, coinctype, ifotime, ifo1, ifo2, ''' H1bgfile ''', nonbgtype, ''' H1nonbgfile ''');'])

  % calculate the ranking for background based on triple det stat
  bgdetstat = sqrt(H1bg.snr.^2 + H2bg.snr.^2 + L1bg.snr.^2);

  % calculate the ranking for pg/zerolag based on triple det stat
  nonbgdetstat = sqrt(H1nonbg.snr.^2 + H2nonbg.snr.^2 + L1nonbg.snr.^2);

  % call another function that plots cumulative hist
  plotcumhist(veto_level, coinctype, ifotime, ifo1, ifo2, bgdetstat, nonbgtype, nonbgdetstat);

  % call another function that plots DT vs F plots for timeslides
  eval(['background_param_plots( veto_level, coinctype, ifotime, ''H1'', ''H2'', ''' H1bgfile ''', ''' H2bgfile ''');'])
  eval(['background_param_plots( veto_level, coinctype, ifotime, ''H1'', ''L1'', ''' H1bgfile ''', ''' L1bgfile ''');'])
  eval(['background_param_plots( veto_level, coinctype, ifotime, ''H2'', ''L1'', ''' H2bgfile ''', ''' L1bgfile ''');'])
