% sampling frequency
dt = 1/16384.0;

% total time
T = 64;

% GPS start time
T0 = 0;

ETG = 'tfclusters';            % ETG to use
channel = 'test';              % channel name
thr = 0.0379269;              % threshold
dowin = 1;                     % window?
thrmethod = 0;                 % threshold method
tfinfo = 0;                    % report tf info?
Tres = 0.03125;                  % time resolution
minf = -130;                      % min frequency
maxf = -400;                   % max frequency
alpha = -1.5;                   % alpha
sigma = 2;                     % sigma
delta = [0,0,0,0,0,0,0,0,0,0,0,0,0]; % delta

% create noise
%N = T ./ dt;
%n = randn(1,N);

[n,t] = frextract('/home/jsylvest/S2Cache/L1-RDS-734338016-16.gwf','L1:LSC-AS_Q');;

% signal
h0 = 1.0;
f0 = 400;
tau = 2/f0;
hp=h0*sin(2*pi*f0*t).*exp(-(t-max(t)/2).^2/tau.^2);

n = n+hp;

% run the ETG
if(strcmp(ETG,'tfclusters'))
 [start_time_s,start_time_ns,duration,central_freq,bandwidth] = LALBurst(n,T0,dt,ETG,channel,'string',thr,'float',dowin,'int',thrmethod,'int',tfinfo,'int',Tres,'float',minf,'float',maxf,'float',alpha,'float',sigma,'int',delta(1),'int',delta(2),'int',delta(3),'int',delta(4),'int',delta(5),'int',delta(6),'int',delta(7),'int',delta(8),'int',delta(9),'int',delta(10),'int');
else
 error('unknown ETG');
end

figure(1)
clf
plot(t,n);
hold on
plot(start_time_s+1e-9*start_time_ns,zeros(size(start_time_s)),'pr');
plot(t,hp,'g')

figure(2)
%clf
hold on
tim = start_time_s + 1e-9*start_time_ns;
for k=1:length(start_time_s)

     h = rectangle('Position',[tim(k) central_freq(k)-0.5*bandwidth(k) duration(k) bandwidth(k)]);
     set(h,'EdgeColor',[0 0 1])

end
