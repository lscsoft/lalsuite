
%
% snr_vs_ds2 is an independent test of the XLAL3DRinca function in
% SnglRingdownUtils.c
%
% It reads in *_trip.mat files that are output from separate.m and
% calculates ds^2_min for each ifo pairing in order to plot SNR vs. ds2
%
% snr_vs_ds2 calls metric3d.m as an external function in order to
% calculate the 3D metric ds^2(f,Q,t)
%
% snr_vs_ds2 then plots and saves SNR(ifo1) vs. ds^2(ifo1,ifo2) for each of
% the 3 detector combinations.
%
% Usage Example:
%
% snr_vs_ds2;
%
% *_trip.mat files and metric3d.m should be in the working directory or
% accesssible path. No other inputs are required.
%
% Author: Fiona Speirits

ifo=['H1';'H2';'L1'];
list=['H1';'H2';'H1';'L1';'H2';'L1']; % all three pairs in order

% read in files
for j=1:3
 eval(['inj' ifo(j,:) 't=load(''inj' ifo(j,:) '_trip.mat'');'])
end

% Loop to select dt resulting in minimum ds^2 to remove signal travel time
% issues for h1l1 and h2l1 -- signal_travel_time_max = 10ms, sampling_rate = 8192Hz

sig_trav_time_max = 1e-2;
sampling_rate = 8192;
step = 1/sampling_rate;

% light travel time 
lightH1H2 = 0;
lightH1L1 = sig_trav_time_max;
lightH2L1 = sig_trav_time_max;

% loop through each of the pairs, calling the function metric3d:
% this follows SnglRingdownUtils very closely
for j=1:2:5
  for i=1:length(injH1t.t)
    eval(['dtab = inj' list(j+1,:) 't.t(i)-inj' list(j,:) 't.t(i);']);
    eval(['dtba = inj' list(j,:) 't.t(i)-inj' list(j+1,:) 't.t(i);']);
    eval(['dt_min_ab = dtab - light' list(j,:) '' list(j+1,:) ';']);
    eval(['dt_min_ba = dtba - light' list(j,:) '' list(j+1,:) ';']);
    eval(['dt_max_ab = dtab + light' list(j,:) '' list(j+1,:) ';']);  
    eval(['dt_max_ba = dtba + light' list(j,:) '' list(j+1,:) ';']);
    dt_min = min(dt_min_ab,dt_min_ba);
    dt_max = max(dt_max_ab,dt_max_ba);
    eval(['ds2_min = metric3d(inj' list(j,:) 't.f(i),inj' list(j+1,:) 't.f(i),inj' list(j,:) 't.Q(i),inj' list(j+1,:) 't.Q(i),dtab);'])
    for dt=dt_min:step:dt_max
      eval(['ds2 = metric3d(inj' list(j,:) 't.f(i),inj' list(j+1,:) 't.f(i),inj' list(j,:) 't.Q(i),inj' list(j+1,:) 't.Q(i),dt);'])
      if ds2 < ds2_min;
        ds2_min = ds2;
      end
    end
    eval(['ds' list(j,:) '' list(j+1,:) '(i)=ds2_min;'])
  end
  % plot snr vs ds^2
  figure
  eval(['semilogy(ds' list(j,:) '' list(j+1,:) ',inj' list(j,:) 't.snr,''.'')'])
  grid on
  eval(['xlabel(''ds^2_{' list(j,:) '' list(j+1,:) '}'')'])
  eval(['ylabel(''SNR_{' list(j,:) '}'')'])
  eval(['saveas(gcf,''snr' list(j,:) 'vds2' list(j,:) '' list(j+1,:) '.png'')'])
end
