%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% triple coincidence data; columns are: time delay(dt), start_time, end_time, central_freq, bandwidth
TripleDat = []; 

% vector of structures for coincidence data
CList = []; 

% counter for triple coincidences
kk=1; 

% path to ETG data that has been processed by the sblit c program; need 3 directories: H1, H2, L1
%path = '/home/jsylvest/S2/S2v2.2/Sblit'; % 130-400 Hz band
%path = '/home/jsylvest/S2/S2.B2.Full/Sblit'; % 400-1000 Hz band

path = '/home/jsylvest/S2/GA2/Sblit/Scratch.-130,-400/';
%path = '/home/jsylvest/S2/GA2/Sblit/Scratch.-400,-1000/';

% Cuts on the SNR for 130-400 Hz band
SNRL1 = 10.6138;
SNRH1 = 3.3113;
SNRH2 = 2.8218;

% Cuts for 400-1000 Hz
%SNRL1 = 5.9738;
%SNRH1 = 5.9738;
%SNRH2 = 5.9738;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try % first check if this script has been run

 s = load2([path '/jobs.mat']); % load mat file

 % data structures
 H1 = s.H1;
 H2 = s.H2;
 L1 = s.L1;

 % time segment structures
 H1Segs = s.H1Segs;
 H2Segs = s.H2Segs;
 L1Segs = s.L1Segs;

catch % do this if first time

 % load all files
 if(exist([path '/H1'],'dir'))
  [H1,H1Segs] = read_burstdso_files([path '/H1/jobH1*.0']);
  [H2,H2Segs] = read_burstdso_files([path '/H2/jobH2*.0']);
  [L1,L1Segs] = read_burstdso_files([path '/L1/jobL1*.0']);
 else
  [H1,H1Segs] = read_burstdso_files([path '/jobH1*.0']);
  [H2,H2Segs] = read_burstdso_files([path '/jobH2*.0']);
  [L1,L1Segs] = read_burstdso_files([path '/jobL1*.0']);
 end

 % generate mat file for future runs
 save([path '/jobs.mat'],'H1','H2','L1','H1Segs','H2Segs','L1Segs');

end % try block


 % peak time: add nanoseconds in
 TL1 = L1.start_time_s+1E-9*L1.start_time_ns;
 TH1 = H1.start_time_s+1E-9*H1.start_time_ns;
 TH2 = H2.start_time_s+1E-9*H2.start_time_ns;

 % sort by peak time
 [TL1,I] = sort(TL1);
 L1d = L1.duration(I);
 L1f = L1.central_freq(I);
 L1b = L1.bandwidth(I);
 L1s = L1.snr(I);

 [TH1,I] = sort(TH1);
 H1d = H1.duration(I);
 H1f = H1.central_freq(I);
 H1b = H1.bandwidth(I);
 H1s = H1.snr(I);

 [TH2,I] = sort(TH2);
 H2d = H2.duration(I);
 H2f = H2.central_freq(I);
 H2b = H2.bandwidth(I);
 H2s = H2.snr(I);

 % segment start times: put nanoseconds in
 L1S = L1Segs.time_s + 1e-9*L1Segs.time_ns;
 H1S = H1Segs.time_s + 1e-9*H1Segs.time_ns;
 H2S = H2Segs.time_s + 1e-9*H2Segs.time_ns;

 % impose SNR cuts
 I = find(L1s >= SNRL1);
 TL1=TL1(I);
 L1d = L1d(I);
 L1f = L1f(I);
 L1s = L1s(I);
 L1b = L1b(I);

 I = find(H1s >= SNRH1);
 TH1=TH1(I);
 H1d = H1d(I);
 H1f = H1f(I);
 H1s = H1s(I);
 H1b = H1b(I);

 I = find(H2s >= SNRH2);
 TH2=TH2(I);
 H2d = H2d(I);
 H2f = H2f(I);
 H2s = H2s(I);
 H2b = H2b(I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coincidence code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NClust = []; % number of coincidences per time shift

dt = -115:5:115; % time shifts, in seconds

Olap = zeros(size(dt)); % uptime per time shift, seconds


% Apply Hanford time-frequency coincidence outside of time shift loop
[IL3,IH3]=timefreq_coin_fancy(TH1-H1d,2*H1d,TH2-H2d,2*H2d,H1f,H2f,H1b,H2b);


for j=1:length(dt) % loop over time shifts
 
  % L1-H1 time coincidence, delay added to H1 times
  % use only H1 events in coincidence with L1
  [IL1,IH1]=time_coin(TL1-L1d,2*L1d,dt(j)+TH1(IL3)-H1d(IL3),2*H1d(IL3));

  % reexpress IH1 in terms of full list
  IH1 = IL3(IH1); 

  % L1-H2 time coincidence, delay added to H2 times
  % use only H2 events in coincidence with L1
  [IL2,IH2]=time_coin(TL1-L1d,2*L1d,dt(j)+TH2(IH3)-H2d(IH3),2*H2d(IH3));

  % reexpress IH2 in terms of full list
  IH2 = IH3(IH2);

  IL1 = intersect(IL1,IL2); % only keep L1 events coincident with H1 and H2

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % form clusters from events that are
  % time coincident
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  tclust = {}; % clusters data
  nclust = 0;  % number of clusters

  % registers, ==1 if event not yet analyzed
  gL1 = ones(size(IL1));
  gH1 = ones(size(IH1));
  gH2 = ones(size(IH2));

  jL1 = 1; % index in L1 events

  while(sum(gL1) > 0) % loop until all L1 events are clustered

     if(jL1 > length(IL1))
      break
     end

     ev.L1 = IL1(jL1); % select a L1 event
     t = TL1(ev.L1(end)); % peak time of event
     du = L1d(ev.L1(end)); % duration of event

     gL1(jL1) = 0; % remember event has been processed

     % containers for H1 and H2 events
     ev.H1 = [];
     ev.H2 = [];

     while(1)
 
      % find Hanford events that are time coincident with L1 event
      i1 = find(gH1 & abs(dt(j)+TH1(IH1)-t) <= du + H1d(IH1));
      i2 = find(gH2 & abs(dt(j)+TH2(IH2)-t) <= du + H2d(IH2));

      % store events
      ev.H1 = [ev.H1;IH1(i1)];
      ev.H2 = [ev.H2;IH2(i2)];

      % test loop
      if(isempty(ev.H1) & isempty(ev.H2))
        jL1 = jL1 + 1;
        break
      end

      % remember events have been processed
      gH1(i1)=0;
      gH2(i2)=0;

      % move to next L1 event
      jL1 = jL1 + 1;

      if(jL1 > length(IL1))
         break
      end

      if(isempty(find(abs(TL1(IL1(jL1)) - TH1(ev.H1) - dt(j)) <= du + H1d(ev.H1)))~=1 | isempty(find(abs(TL1(IL1(jL1)) - TH2(ev.H2) - dt(j)) <= du + H2d(ev.H2)))~=1)
        % if L1 event coincident with any of H events in cluster, store and continue
        ev.L1 = [ev.L1 ; IL1(jL1)];
        t = TL1(ev.L1(end));
	du = L1d(ev.L1(end));
        gL1(jL1) = 0;
      else
	% otherwise, break
        break
      end

     end % while loop


     nclust = nclust + 1; % increase number of clusters
     tclust{nclust} = ev; % save cluster data

end % while(sum(gL1) > 0)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge clusters within ctime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     ctime = 0.2; % clustering timescale, seconds

     for k=2:length(tclust) % loop over clusters

      if(min([TL1(tclust{k}.L1);dt(j)+TH1(tclust{k}.H1);dt(j)+TH2(tclust{k}.H2)]) - max([TL1(tclust{k-1}.L1);dt(j)+TH1(tclust{k-1}.H1);dt(j)+TH2(tclust{k-1}.H2)]) <= ctime)
 
       % clusters overlap; merge
       tclust{k}.L1 = [tclust{k-1}.L1 ; tclust{k}.L1];
       tclust{k}.H1 = [tclust{k-1}.H1 ; tclust{k}.H1];
       tclust{k}.H2 = [tclust{k-1}.H2 ; tclust{k}.H2];
       tclust{k-1}.L1 = [];
       tclust{k-1}.H1 = [];
       tclust{k-1}.H2 = [];

      end % if
 
     end % for

     % expunge empty clusters
     nclust = 0;
     tc2 = {};
     for k=1:length(tclust)
      if(isempty(tclust{k}.L1) ~= 1)
       nclust = nclust+1;
       tc2{nclust} = tclust{k};
      end
     end
     tclust = tc2;

     % impose strict coincidence within clusters
     good=[]; % container 

     for k=1:length(tclust) % loop over clusters

        ev = tclust{k}; % a cluster

        good(k) = 0; % cluster is assumed bad

        for n=1:length(ev.L1) % loop over L1 events in cluster

	 % find H events coincident with L1 event
         I1 = find(abs(TL1(ev.L1(n))-TH1(ev.H1)-dt(j)) <= L1d(ev.L1(n)) + H1d(ev.H1));
         I2 = find(abs(TL1(ev.L1(n))-TH2(ev.H2)-dt(j)) <= L1d(ev.L1(n)) + H2d(ev.H2));

	 for m=1:length(I1) % check consistency of H events
          I3 = find(abs(TH1(ev.H1(I1(m))) - TH2(ev.H2(I2))) <= H1d(ev.H1(I1(m))) + H2d(ev.H2(I2)));

          if(length(I3)>0)
            % have at least one set of triple coincident events in cluster
            good(k) = 1;
          end

         end % for m=1:length(I1)
       end % for n=1:length(ev.L1)
     end % for k=1:length(tclust)


     % Save number of clusters for time shift j
     NClust(j) = sum(good);

     disp(num2str([dt(j) NClust(j)]));


     %%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Calculate overlap
     %%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     % find segment times HT with H1 and H2 data
     HT = [];
     for k=1:length(H1S)
      if(isempty(find(H1S(k) == H2S))~=1) 
       HT=[HT;H1S(k)];
      end
     end

     for k=1:length(L1S) % loop over L1 segment times

      % times overlaping with HT
      I = find(L1S(k) + 300 >= HT + dt(j) & L1S(k) <= HT + 300 + dt(j));

      for m=1:length(I) % loop over overlaping times

       % Calculate contribution to overlap
       if(HT(I(m))+dt(j) >= L1S(k))
        Olap(j) = Olap(j) + L1S(k) + 300 - HT(I(m))-dt(j);
       else
        Olap(j) = Olap(j) + HT(I(m)) + 300 - L1S(k)+dt(j);
       end

      end % for m=1:length(I)

     end % for k=1:length(L1S)


     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Build final event list
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     for k=1:length(tclust) % loop over clusters

      if(good(k)) % if they are good:

       % start time = mean of extrema
       CList(kk).time = (min([TL1(tclust{k}.L1);dt(j)+TH1(tclust{k}.H1);dt(j)+TH2(tclust{k}.H2)]) + max([TL1(tclust{k}.L1);dt(j)+TH1(tclust{k}.H1);dt(j)+TH2(tclust{k}.H2)])) / 2;

       % duration = difference of extrema
       CList(kk).duration = max([TL1(tclust{k}.L1);dt(j)+TH1(tclust{k}.H1);dt(j)+TH2(tclust{k}.H2)]) - min([TL1(tclust{k}.L1);dt(j)+TH1(tclust{k}.H1);dt(j)+TH2(tclust{k}.H2)]);

       % central_freq = mean of extrema
       CList(kk).central_freq = (min([L1f(tclust{k}.L1);H1f(tclust{k}.H1);H2f(tclust{k}.H2)]) + max([L1f(tclust{k}.L1);H1f(tclust{k}.H1);H2f(tclust{k}.H2)])) / 2;

       % bandwidth = difference of extrema
       CList(kk).bandwidth = max([L1f(tclust{k}.L1);H1f(tclust{k}.H1);H2f(tclust{k}.H2)]) - min([L1f(tclust{k}.L1);H1f(tclust{k}.H1);H2f(tclust{k}.H2)]);

       disp(num2str([dt(j) CList(kk).time CList(kk).time+CList(kk).duration CList(kk).central_freq CList(kk).bandwidth]));

       % store for rStatistics
       TripleDat = [TripleDat ; [dt(j) CList(kk).time CList(kk).time+CList(kk).duration CList(kk).central_freq CList(kk).bandwidth]];

    
       kk = kk + 1;


       end % if good
       
     end % loop over clusters

end % loop over time shifts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% events
Dat = [TripleDat(:,2) TripleDat(:,3)-TripleDat(:,2) TripleDat(:,1)];

fid = fopen([path '/triple_coin_S2.dat'],'w');

fprintf(fid,'%.15g\t%g\t0\t0\t%g\n',Dat');

fclose(fid);

% uptime
fid = fopen([path '/livetime.dat'],'w');

fprintf(fid,'%g\t%g\n',[dt' Olap']');

fclose(fid)
