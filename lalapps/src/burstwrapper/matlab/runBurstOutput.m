Ntrial = 1;

%noise
N = 64 * 16384;
t0 = 0;
dt = 1/16384;

%signal
A = 0.1;
t = 0:1/16384:1-1/16384;
cf = 250;
Q = 9;


% define response function
response_df = 0.1; % freq resolution
response_f0 = 0; % starting frequency
response_re = ones(1,ceil(0.5/dt/response_df));
response_im = zeros(1,ceil(0.5/dt/response_df));

tbias = [];
tvar = [];
fbias = [];
fvar = [];
msnr = [];
mamp = [];
vsnr = [];
vamp = [];

for jj=1:length(A)

     disp(['Now doing: A = ' num2str(A(jj))]);

% run
terr = [];
ferr = [];
snr = [];
amp = [];
con = [];

for j=1:Ntrial

if(mod(j,floor(Ntrial/10))==0)
     disp(['.........' num2str(100*j/Ntrial,2) '%']);
end

n = randn(1,N);
dts = Q/(sqrt(2)*pi*cf);
sig = A(jj)*sin(2*pi*cf*(t-0.5)).*exp(-(t-0.5).^2/(dts^2))/sqrt(pi*dts^2);
sigsnr = sqrt(sum(sig.^2));

sig0 = [ zeros(1,(N-length(sig))/2) sig zeros(1,(N-length(sig))/2) ];
y = n + sig0;

duration = 4/16;
start_time_s = floor(t0 + dt * N / 2 - duration / 2);
start_time_ns = 1E9*(t0 + dt * N / 2 - duration / 2 - start_time_s);
central_freq = cf;
bandwidth = 200;

[junk,maxtime] = max(sig0);
maxtime = t0 + dt * maxtime;

[start_time_s_,start_time_ns_,duration_,central_freq_,bandwidth_,amplitude_,snr_,confidence_] = LALBurstOutput(y,t0,dt,response_re,response_im,response_f0,response_df,start_time_s,start_time_ns,duration,central_freq,bandwidth);

terr = [ terr ; maxtime - start_time_s_ - 1e-9*start_time_ns_ ];
ferr = [ ferr ; central_freq - central_freq_ ];
snr = [ snr ; snr_ ];
amp = [ amp ; amplitude_ ];
con = [ con ; confidence_ ];


plot(sig0);
hold on;
plot((start_time_s_+1E-9*start_time_ns_)/dt,0,'r.');

drawnow

end

tbias(jj) = mean(terr);
tvar(jj) = mean((terr-tbias(jj)).^2);

fbias(jj) = mean(ferr);
fvar(jj) = mean((ferr-fbias(jj)).^2);

msnr(jj) = mean(snr);
vsnr(jj) = var(snr);

mamp(jj) = mean(amp);
vamp(jj) = var(amp);

end

%save('/home/jsylvest/c/testBurstOutput.mat','A','tbias','tvar','fbias','fvar','msnr','vsnr','mamp','vamp');
