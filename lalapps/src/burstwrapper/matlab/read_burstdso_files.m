function [Dat,Segs] = read_burstdso_files(path,I,username);

if(nargin<=2)
     [st,username] = system('whoami');
     if(st) 
     error('cannot get user name from whoami; use read_burstdso_files(path,I,username)');
     end
end

if(isempty(strfind(path,'*')))
 dd = dir([path '/job*']);
else
 dd = dir(path);
 ii = strfind(path,'/');
 ii = ii(end);
 path = path(1:ii);
end

Dat=struct('start_time_s',[],'start_time_ns',[],'duration',[],'central_freq',[],'bandwidth',[],'amplitude',[],'snr',[],'confidence',[]);

Segs = struct('time_s',[],'time_ns',[],'Nevents',[]);

if(nargin == 1)
     I = 1:length(dd);
end

for jj=1:length(I)
     j = I(jj);
     disp(dd(j).name);

     scanline = ['job%c%i.' username '.%i.%i.bin.%i.%i.0'];
     A = sscanf(dd(j).name,scanline);
     Segs.time_s = [Segs.time_s ; A(5)];
     Segs.time_ns = [Segs.time_ns ; A(6)];

     [ifo,search,channel,tstart_time_s,tstart_time_ns,tduration,tcentral_freq,tbandwidth,tamplitude,tsnr,tconfidence] = reab([path '/' dd(j).name]);
     Dat.start_time_s = [Dat.start_time_s ; tstart_time_s];
     Dat.start_time_ns = [Dat.start_time_ns ; tstart_time_ns];
     Dat.duration = [Dat.duration;tduration];
     Dat.central_freq = [Dat.central_freq;tcentral_freq];
     Dat.bandwidth = [Dat.bandwidth;tbandwidth];
     Dat.amplitude = [Dat.amplitude;tamplitude];
     Dat.snr = [ Dat.snr;tsnr ];
     Dat.confidence = [Dat.confidence;tconfidence];

     Segs.Nevents = [ Segs.Nevents ; length(tsnr) ];
end
