% reae are very few injections, but why do these look more symmetric
%clear

% read in timeslide arrays
% cd timeslides
%{
bgH1t=load('bgH1_trip.mat');
bgL1t=load('bgL1_trip.mat');
bgH2t=load('bgH2_trip.mat');
bgH1inL1d=load('bgH1inL1_doub.mat');
bgL1inH1d=load('bgL1inH1_doub.mat');
bgH1inH2d=load('bgH1inH2_doub.mat');
bgH2inH1d=load('bgH2inH1_doub.mat');
bgL1inH2d=load('bgL1inH2_doub.mat');
bgH2inL1d=load('bgH2inL1_doub.mat');
%}

% read in the injection lists
injH1t=load('injH1_trip.mat');
injL1t=load('injL1_trip.mat');
injH2t=load('injH2_trip.mat');
injH1d=load('injH1d.mat');
injH2d=load('injH2d.mat');
injL1d=load('injL1d.mat');
injH1inL1d=load('injH1inL1_doub.mat');
injL1inH1d=load('injL1inH1_doub.mat');
injH1inH2d=load('injH1inH2_doub.mat');
injH2inH1d=load('injH2inH1_doub.mat');
injL1inH2d=load('injL1inH2_doub.mat');
injH2inL1d=load('injH2inL1_doub.mat');


%=========================== Raise doubles threshold ========================%

% cut on loudest event stat 

%from background estimate
%{
injH1inL1dcut.ind = injH1inL1d.ind(injH1inL1d.dst<20.4);
injH1inL1daftercut.ind = setdiff(injH1inL1d.ind,injH1inL1dcut.ind);
injL1inH1dcut.ind = injL1inH1d.ind(injL1inH1d.dst<20.4);
injL1inH1daftercut.ind = setdiff(injL1inH1d.ind,injL1inH1dcut.ind);

injH1inH2dcut.ind = injH1inH2d.ind(injH1inH2d.dst<15.07);
injH1inH2daftercut.ind = setdiff(injH1inH2d.ind,injH1inH2dcut.ind);
injH2inH1dcut.ind = injH2inH1d.ind(injH2inH1d.dst<15.07);
injH2inH1daftercut.ind = setdiff(injH2inH1d.ind,injH2inH1dcut.ind);

injH2inL1dcut.ind = injH2inL1d.ind(injH2inL1d.dst<18.85);
injH2inL1daftercut.ind = setdiff(injH2inL1d.ind,injH2inL1dcut.ind);
injL1inH2dcut.ind = injL1inH2d.ind(injL1inH2d.dst<18.85);
injL1inH2daftercut.ind = setdiff(injL1inH2d.ind,injL1inH2dcut.ind);
%}

% loudest intime event
%{
injH1inL1dcut.ind = injH1inL1d.ind(injH1inL1d.dst<21.05);
injH1inL1daftercut.ind = setdiff(injH1inL1d.ind,injH1inL1dcut.ind);
injL1inH1dcut.ind = injL1inH1d.ind(injL1inH1d.dst<21.05); 
injL1inH1daftercut.ind = setdiff(injL1inH1d.ind,injL1inH1dcut.ind); 

injH1inH2dcut.ind = injH1inH2d.ind(injH1inH2d.dst<63.10); 
injH1inH2daftercut.ind = setdiff(injH1inH2d.ind,injH1inH2dcut.ind); 
injH2inH1dcut.ind = injH2inH1d.ind(injH2inH1d.dst<63.10);
injH2inH1daftercut.ind = setdiff(injH2inH1d.ind,injH2inH1dcut.ind);

injH2inL1dcut.ind = injH2inL1d.ind(injH2inL1d.dst<20.27);
injH2inL1daftercut.ind = setdiff(injH2inL1d.ind,injH2inL1dcut.ind);
injL1inH2dcut.ind = injL1inH2d.ind(injL1inH2d.dst<20.27);
injL1inH2daftercut.ind = setdiff(injL1inH2d.ind,injL1inH2dcut.ind);
%}


% if no cut

injH1inL1dcut.ind = [];
injH1inH2dcut.ind = [];
injH2inL1dcut.ind = [];
injH1inL1daftercut.ind = injH1inL1d.ind;
injH1inH2daftercut.ind = injH1inH2d.ind;
injH2inL1daftercut.ind = injH2inL1d.ind;


%=========================== Injected Quantities ===========================%

% Read in missed and found sim_ringdown tables

%A=[2,3,4,5,7,8,9,10,11];
A=1;
eval(['fsim=readMeta(''injH1H2L1coincs.xml'',''sim_ringdown'',0,''h_start_time,h_start_time_ns,l_start_time,l_start_time_ns,mass,spin,frequency,quality,eff_dist_h,eff_dist_l,distance,hrss_h,hrss_l'');'])

eval(['msim=readMeta(''injH1H2L1missed.xml'',''sim_ringdown'',0,''h_start_time,h_start_time_ns,l_start_time,l_start_time_ns,mass,spin,frequency,quality,eff_dist_h,eff_dist_l,distance,hrss_h,hrss_l'');'])

for k=1:length(fsim.mass)
  fsim.run(k)=A(1);
end
fsim.run=transpose(fsim.run);

for k=1:length(msim.mass)
  msim.run(k)=A(1);
end
msim.run=transpose(msim.run);

%{
for i=2:9
eval(['fsimi=readMeta(''H1H2L1coincs_' num2str(A(i)) '-h1h2cl.xml'',''sim_ringdown'',0,''h_start_time,h_start_time_ns,l_start_time,l_start_time_ns,mass,spin,quality,frequency,eff_dist_h,eff_dist_l,distance,hrss_h,hrss_l'');'])
eval(['msimi=readMeta(''H1H2L1missed_' num2str(A(i)) '-h1h2cl.xml'',''sim_ringdown'',0,''h_start_time,h_start_time_ns,l_start_time,l_start_time_ns,mass,spin,quality,frequency,eff_dist_h,eff_dist_l,distance,hrss_h,hrss_l'');'])

  for k=1:length(fsimi.quality)
    fsimi.run(k)=A(i);
  end
  fsimi.run=transpose(fsimi.run);

  for k=1:length(msimi.quality)
    msimi.run(k)=A(i);
  end
  msimi.run=transpose(msimi.run);

  fsim.h_start_time=[fsim.h_start_time;fsimi.h_start_time];
  fsim.h_start_time_ns=[fsim.h_start_time_ns;fsimi.h_start_time_ns];
  fsim.l_start_time=[fsim.l_start_time;fsimi.l_start_time];
  fsim.l_start_time_ns=[fsim.l_start_time_ns;fsimi.l_start_time_ns];
  fsim.frequency=[fsim.frequency;fsimi.frequency];
  fsim.quality=[fsim.quality;fsimi.quality];
  fsim.mass=[fsim.mass;fsimi.mass];  
  fsim.spin=[fsim.spin;fsimi.spin];
  fsim.eff_dist_h=[fsim.eff_dist_h;fsimi.eff_dist_h];
  fsim.eff_dist_l=[fsim.eff_dist_l;fsimi.eff_dist_l];
  fsim.distance=[fsim.distance;fsimi.distance];
  fsim.run=[fsim.run;fsimi.run];
  fsim.hrss_h=[fsim.hrss_h;fsimi.hrss_h];
  fsim.hrss_l=[fsim.hrss_l;fsimi.hrss_l];
  msim.h_start_time=[msim.h_start_time;msimi.h_start_time];
  msim.h_start_time_ns=[msim.h_start_time_ns;msimi.h_start_time_ns];
  msim.l_start_time=[msim.l_start_time;msimi.l_start_time];
  msim.l_start_time_ns=[msim.l_start_time_ns;msimi.l_start_time_ns];
  msim.frequency=[msim.frequency;msimi.frequency];
  msim.quality=[msim.quality;msimi.quality];
  msim.mass=[msim.mass;msimi.mass];           
  msim.spin=[msim.spin;msimi.spin];
  msim.eff_dist_h=[msim.eff_dist_h;msimi.eff_dist_h];
  msim.eff_dist_l=[msim.eff_dist_l;msimi.eff_dist_l];
  msim.distance=[msim.distance;msimi.distance];
  msim.run=[msim.run;msimi.run];
  msim.hrss_h=[msim.hrss_h;msimi.hrss_h];
  msim.hrss_l=[msim.hrss_l;msimi.hrss_l];

end
%}
msim.th=msim.h_start_time+msim.h_start_time_ns/1e9;
msim.tl=msim.l_start_time+msim.l_start_time_ns/1e9;
msim.f=msim.frequency;
msim.q=msim.quality;
msim.m=msim.mass;
msim.a=msim.spin;
msim.d=msim.distance;
msim.dh=msim.eff_dist_h;
msim.dl=msim.eff_dist_l;

fsim.th=fsim.h_start_time+fsim.h_start_time_ns/1e9;
fsim.tl=fsim.l_start_time+fsim.l_start_time_ns/1e9;
fsim.f=fsim.frequency;
fsim.q=fsim.quality;
fsim.m=fsim.mass;
fsim.a=fsim.spin;
fsim.d=fsim.distance;
fsim.dh=fsim.eff_dist_h;
fsim.dl=fsim.eff_dist_l;

for i=1:length(fsim.th)
  fsim.ind(i)=i;
end

% the index of fsim is the same as the fsim.ind
fsimH1L1d.ind=fsim.ind(injH1inL1daftercut.ind);
fsimH1L1d.f=fsim.f(injH1inL1daftercut.ind);
fsimH1L1d.d=fsim.d(injH1inL1daftercut.ind);
fsimH1L1d.dh=fsim.dh(injH1inL1daftercut.ind);
fsimH1L1d.dl=fsim.dl(injH1inL1daftercut.ind);
fsimH1L1d.th=fsim.th(injH1inL1daftercut.ind);
fsimH1L1d.tl=fsim.tl(injH1inL1daftercut.ind);
fsimH1L1d.run=fsim.run(injH1inL1daftercut.ind);
fsimH1L1d.hrss_h=fsim.hrss_h(injH1inL1daftercut.ind);
fsimH1L1d.hrss_l=fsim.hrss_l(injH1inL1daftercut.ind);

fsimH1H2d.ind=fsim.ind(injH1inH2daftercut.ind);
fsimH1H2d.f=fsim.f(injH1inH2daftercut.ind);
fsimH1H2d.d=fsim.d(injH1inH2daftercut.ind);
fsimH1H2d.dh=fsim.dh(injH1inH2daftercut.ind);
fsimH1H2d.dl=fsim.dl(injH1inH2daftercut.ind);
fsimH1H2d.th=fsim.th(injH1inH2daftercut.ind);
fsimH1H2d.tl=fsim.tl(injH1inH2daftercut.ind);
fsimH1H2d.run=fsim.run(injH1inH2daftercut.ind);
fsimH1H2d.hrss_h=fsim.hrss_h(injH1inH2daftercut.ind);
fsimH1H2d.hrss_l=fsim.hrss_l(injH1inH2daftercut.ind);

fsimL1H2d.ind=fsim.ind(injH2inL1daftercut.ind);
fsimL1H2d.f=fsim.f(injH2inL1daftercut.ind);
fsimL1H2d.d=fsim.d(injH2inL1daftercut.ind);
fsimL1H2d.dh=fsim.dh(injH2inL1daftercut.ind);
fsimL1H2d.dl=fsim.dl(injH2inL1daftercut.ind);
fsimL1H2d.th=fsim.th(injH2inL1daftercut.ind);
fsimL1H2d.tl=fsim.tl(injH2inL1daftercut.ind);
fsimL1H2d.run=fsim.run(injH2inL1daftercut.ind);
fsimL1H2d.hrss_h=fsim.hrss_h(injH2inL1daftercut.ind);
fsimL1H2d.hrss_l=fsim.hrss_l(injH2inL1daftercut.ind);

foundd=[fsimH1L1d.ind,fsimH1H2d.ind,fsimL1H2d.ind];
threshcut.ind=fsim.ind([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);
threshcut.f=fsim.f([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);
threshcut.d=fsim.d([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);
threshcut.dh=fsim.dh([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);
threshcut.dl=fsim.dl([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);
threshcut.th=fsim.th([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);
threshcut.tl=fsim.tl([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);
threshcut.run=fsim.run([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);
threshcut.hrss_h=fsim.hrss_h([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);
threshcut.hrss_l=fsim.hrss_l([injH1inH2dcut.ind,injH2inL1dcut.ind,injH1inL1dcut.ind]);

msim.th=[msim.th;threshcut.th];
msim.tl=[msim.tl;threshcut.tl];
msim.f=[msim.f;threshcut.f];
msim.d=[msim.d;threshcut.d];
msim.dh=[msim.dh;threshcut.dh];
msim.dl=[msim.dl;threshcut.dl];
msim.run=[msim.run;threshcut.run];
msim.hrss_h=[msim.hrss_h;threshcut.hrss_h];
msim.hrss_l=[msim.hrss_l;threshcut.hrss_l];

% read in Veto files
%{
%{
H1veto=load('H1vetocat23.txt');
H1vetostart=H1veto(:,1);
H1vetoend=H1veto(:,2);
% make a list of all the gps seconds vetoed
j=1;
for i=1:length(H1vetostart)
  for k=H1vetostart(i):H1vetoend(i)-1
    H1vetolist(j)=k;
    j=j+1;
  end
end
% save as a mat file
save H1vetolist.mat H1vetolist
%}
A=load('H1vetolist.mat');
H1vetolist=A.H1vetolist;

%{
H2veto=load('H2vetocat23.txt');
H2vetostart=H2veto(:,1);
H2vetoend=H2veto(:,2);
j=1;
for i=1:length(H2vetostart)
  for k=H2vetostart(i):H2vetoend(i)-1
    H2vetolist(j)=k;
    j=j+1;
  end
end
save H2vetolist.mat H2vetolist
%}
B=load('H2vetolist.mat');
H2vetolist=B.H2vetolist;

%{
L1veto=load('L1vetocat23.txt');
L1vetostart=L1veto(:,1);
L1vetoend=L1veto(:,2);
j=1;
for i=1:length(L1vetostart)
  for k=L1vetostart(i):L1vetoend(i)-1
    L1vetolist(j)=k;
    j=j+1;
  end
end
save L1vetolist.mat L1vetolist
%}
C=load('L1vetolist.mat');
L1vetolist=C.L1vetolist;
%}
H1vetolist=[];
H2vetolist=[];
L1vetolist=[];

% times vetoed in all three ifos
vetoedinH1H2L1=intersect(H1vetolist,intersect(H2vetolist,L1vetolist));
% times vetoed in two ifos (only)
vetoedinH1H2=setdiff(intersect(H1vetolist,H2vetolist),vetoedinH1H2L1);
vetoedinH1L1=setdiff(intersect(H1vetolist,L1vetolist),vetoedinH1H2L1);
vetoedinL1H2=setdiff(intersect(L1vetolist,H2vetolist),vetoedinH1H2L1);
% time vetoed in one ifo (only)
vetoedinH1=setdiff(H1vetolist,union(H2vetolist,L1vetolist));
vetoedinH2=setdiff(H2vetolist,union(H1vetolist,L1vetolist));
vetoedinL1=setdiff(L1vetolist,union(H2vetolist,H1vetolist));

% total triple time lost
losttriptimes=[vetoedinH1H2L1,vetoedinH1H2,vetoedinH1L1,vetoedinL1H2,vetoedinH1,vetoedinH2,vetoedinL1];

% list the injections that were missed during these times
vetot.t=msim.th(ismember(floor(msim.th),losttriptimes));
% list the injections that were missed outside of these times (ie MISSED in TRIPLE time)
[missed.th,ind]=setdiff(msim.th,vetot.t);
missed.f=msim.f(ind);
missed.d=msim.d(ind);
missed.dh=msim.dh(ind);
missed.dl=msim.dl(ind);
missed.run=msim.run(ind);
missed.hrss_h=msim.hrss_h(ind);
missed.hrss_l=msim.hrss_l(ind);
clear ind 


% list the injections that were FOUND (as DOUBLES) in DOUBLE time

fsimL1H2dd.th=fsimL1H2d.th(ismember(floor(fsimL1H2d.th),vetoedinH1));
fsimL1H2dd.tl=fsimL1H2d.tl(ismember(floor(fsimL1H2d.th),vetoedinH1));
fsimL1H2dd.d=fsimL1H2d.d(ismember(floor(fsimL1H2d.th),vetoedinH1));
fsimL1H2dd.dh=fsimL1H2d.dh(ismember(floor(fsimL1H2d.th),vetoedinH1));
fsimL1H2dd.dl=fsimL1H2d.dl(ismember(floor(fsimL1H2d.th),vetoedinH1));
fsimL1H2dd.f=fsimL1H2d.f(ismember(floor(fsimL1H2d.th),vetoedinH1));
fsimL1H2dd.run=fsimL1H2d.run(ismember(floor(fsimL1H2d.th),vetoedinH1));
fsimL1H2dd.hrss_h=fsimL1H2d.hrss_h(ismember(floor(fsimL1H2d.th),vetoedinH1));
fsimL1H2dd.hrss_l=fsimL1H2d.hrss_l(ismember(floor(fsimL1H2d.th),vetoedinH1));

fsimH1L1dd.th=fsimH1L1d.th(ismember(floor(fsimH1L1d.th),vetoedinH2));
fsimH1L1dd.tl=fsimH1L1d.tl(ismember(floor(fsimH1L1d.th),vetoedinH2));
fsimH1L1dd.d=fsimH1L1d.d(ismember(floor(fsimH1L1d.th),vetoedinH2));
fsimH1L1dd.dh=fsimH1L1d.dh(ismember(floor(fsimH1L1d.th),vetoedinH2));
fsimH1L1dd.dl=fsimH1L1d.dl(ismember(floor(fsimH1L1d.th),vetoedinH2));
fsimH1L1dd.f=fsimH1L1d.f(ismember(floor(fsimH1L1d.th),vetoedinH2));
fsimH1L1dd.run=fsimH1L1d.run(ismember(floor(fsimH1L1d.th),vetoedinH2));
fsimH1L1dd.hrss_h=fsimH1L1d.hrss_h(ismember(floor(fsimH1L1d.th),vetoedinH2));
fsimH1L1dd.hrss_l=fsimH1L1d.hrss_l(ismember(floor(fsimH1L1d.th),vetoedinH2));

fsimH1H2dd.th=fsimH1H2d.th(ismember(floor(fsimH1H2d.th),vetoedinL1));
fsimH1H2dd.tl=fsimH1H2d.tl(ismember(floor(fsimH1H2d.th),vetoedinL1));
fsimH1H2dd.d=fsimH1H2d.d(ismember(floor(fsimH1H2d.th),vetoedinL1));
fsimH1H2dd.dh=fsimH1H2d.dh(ismember(floor(fsimH1H2d.th),vetoedinL1));
fsimH1H2dd.dl=fsimH1H2d.dl(ismember(floor(fsimH1H2d.th),vetoedinL1));
fsimH1H2dd.f=fsimH1H2d.f(ismember(floor(fsimH1H2d.th),vetoedinL1));
fsimH1H2dd.run=fsimH1H2d.run(ismember(floor(fsimH1H2d.th),vetoedinL1));
fsimH1H2dd.hrss_h=fsimH1H2d.hrss_h(ismember(floor(fsimH1H2d.th),vetoedinL1));
fsimH1H2dd.hrss_l=fsimH1H2d.hrss_l(ismember(floor(fsimH1H2d.th),vetoedinL1));

fdoubs=[fsimL1H2dd.th;fsimH1L1dd.th;fsimH1H2dd.th];
doubtimes=[vetoedinH1,vetoedinH2,vetoedinL1];

% list the injections that were  MISSED in DOUBLE time
missedd.th=msim.th(ismember(floor(msim.th),doubtimes));
missedd.tl=msim.tl(ismember(floor(msim.th),doubtimes));
missedd.f=msim.f(ismember(floor(msim.th),doubtimes));
missedd.d=msim.d(ismember(floor(msim.th),doubtimes));
missedd.dh=msim.dh(ismember(floor(msim.th),doubtimes));
missedd.dl=msim.dl(ismember(floor(msim.th),doubtimes));
missedd.run=msim.run(ismember(floor(msim.th),doubtimes));
missedd.hrss_h=msim.hrss_h(ismember(floor(msim.th),doubtimes));
missedd.hrss_l=msim.hrss_l(ismember(floor(msim.th),doubtimes));

% list the injections that were FOUND as DOUBLES in TRIPLE time only
[num,ind]=setdiff(fsimH1L1d.th,fsimH1L1dd.th);
fsimH1L1dt.th=fsimH1L1d.th(ind);
fsimH1L1dt.tl=fsimH1L1d.tl(ind);
fsimH1L1dt.f=fsimH1L1d.f(ind);
fsimH1L1dt.d=fsimH1L1d.d(ind);
fsimH1L1dt.dh=fsimH1L1d.dh(ind);
fsimH1L1dt.dl=fsimH1L1d.dl(ind);
fsimH1L1dt.run=fsimH1L1d.run(ind);
fsimH1L1dt.hrss_h=fsimH1L1d.hrss_h(ind);
fsimH1L1dt.hrss_l=fsimH1L1d.hrss_l(ind);
clear num ind

[num,ind]=setdiff(fsimH1H2d.th,fsimH1H2dd.th);
fsimH1H2dt.th=fsimH1H2d.th(ind);
fsimH1H2dt.tl=fsimH1H2d.tl(ind);
fsimH1H2dt.f=fsimH1H2d.f(ind);
fsimH1H2dt.d=fsimH1H2d.d(ind);
fsimH1H2dt.dh=fsimH1H2d.dh(ind);
fsimH1H2dt.dl=fsimH1H2d.dl(ind);
fsimH1H2dt.run=fsimH1H2d.run(ind);
fsimH1H2dt.hrss_h=fsimH1H2d.hrss_h(ind);
fsimH1H2dt.hrss_l=fsimH1H2d.hrss_l(ind);
clear num ind

[num,ind]=setdiff(fsimL1H2d.th,fsimL1H2dd.th);
fsimL1H2dt.th=fsimL1H2d.th(ind);
fsimL1H2dt.tl=fsimL1H2d.tl(ind);
fsimL1H2dt.f=fsimL1H2d.f(ind);
fsimL1H2dt.d=fsimL1H2d.d(ind);
fsimL1H2dt.dh=fsimL1H2d.dh(ind);
fsimL1H2dt.dl=fsimL1H2d.dl(ind);
fsimL1H2dt.run=fsimL1H2d.run(ind);
fsimL1H2dt.hrss_h=fsimL1H2d.hrss_h(ind);
fsimL1H2dt.hrss_l=fsimL1H2d.hrss_l(ind);

%totals
doubsintrip.th=[fsimH1L1dt.th;fsimH1H2dt.th;fsimL1H2dt.th];

% Bookkeeping
% The number of missed and found injections in each category should equal the total number 
% of missed and found injections read in
total=length(msim.th(ismember(floor(msim.th),vetoedinH1H2L1)))... % missed because of H1H2L1 veto
     +length(msim.th(ismember(floor(msim.th),vetoedinH1H2)))...   % missed because of H1H2 veto
     +length(msim.th(ismember(floor(msim.th),vetoedinH1L1)))...   % missed because of H1L1 veto
     +length(msim.th(ismember(floor(msim.th),vetoedinL1H2)))...   % missed because of H2L1 veto
     +length(msim.th(ismember(floor(msim.th),vetoedinH1)))...     % missed in double time
     +length(msim.th(ismember(floor(msim.th),vetoedinH2)))...     % missed in double time
     +length(msim.th(ismember(floor(msim.th),vetoedinL1)))...     % missed in double time (missedd.th etc)
     +length(missed.f)...                                       % missed in triple time
     +length(injH1t.t)...                                       % found in triple time
     +length(doubsintrip.th)...                                  % found doubles in triple time
     +length(fdoubs);                                         % found doubles in double time

readin=length(fsim.th)-length(injH1inL1dcut.ind )-length(injH2inL1dcut.ind )-length(injH1inH2dcut.ind )+length(msim.th);


