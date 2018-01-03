% making up 0.25 % UL using the old run

%prefix='/local_data/sintes/S4/MultiMC_100_199_500/MCfreq_';
%prefix='/local_data/sintes/S4/MultiMC/MultiMC_100_599_500/MCfreq_';
prefix='/local_data/sintes/S4/MultiMC/MultiMC_600_999.5_500/MCfreq_';

%fileoutput = 'MultiMC_100_600.25_coarse';
fileoutput = 'MultiMC_600_1000.25_coarse';

fid = fopen(fileoutput, 'w');

% file with driver output
file = 'SigMax25.mat';
load(file);

Bands = BandList(:,1);
%Nbands = length(Bands);

%The first 100-600Hz
%Nbands = 2000; %valid for 100-600 Hz

%The 600-1000 Hz
Nbands = 1600; %valid for 600-1000 Hz

%fshift=0;   %valid for 100-600 Hz
fshift=2000; %valid for 600-1000 Hz

for bandnumber = fshift+1:(fshift+Nbands); %the  current frequency band
   fmin = BandList(bandnumber, 1);
   fmax = BandList(bandnumber, 2);
   
   %fminNum= floor(fmin); %valid for 100-600 Hz
   fminNum= floor(fmin*2)/2; %valid for 600-1000 Hz
   
   basestring = strcat(prefix, num2str( fminNum ) );
   h0string = strcat(basestring, '_h0');
   ncstring = strcat(basestring, '_nc');
   Ncount = load(ncstring);
   h0val  = load(h0string);
   
   nMonteCarlos=length(Ncount);
   nh0=length(h0val);
   
   fprintf(fid,'%d  %d %d ', fmin, fmax, Sigmax(bandnumber) );
   
   
   for h0num=1:nh0
      x=Ncount(:, h0num+1);
      kkcount = find(x>Sigmax(bandnumber));
      CH(h0num) = length(kkcount)/nMonteCarlos;
      fprintf(fid,' %d %d ', h0val(h0num), CH(h0num) );
   end
   

   h0vec =h0val;
   CLvec = CH;
   UL=0;

  indices = find(CLvec > 0.95);
  if( length(indices)>0)
    k=indices(1);
    h02=h0vec(k);
    CL2=CLvec(k);
    if(k>1)
      h01=h0vec(k-1);
      CL1=CLvec(k-1);
    else
      h01 = h0vec(1);
      CL1 = CLvec(1);
      h02 = h0vec(2);
      CL2 = CLvec(2);
    end
      slope = (h02-h01)/(CL2-CL1);
      UL=h01 +slope*(0.95 -CL1);
  else
    indices = find(CLvec < 0.95);
    k = length(indices);
    h02 = h0vec(k);
    CL2 = CLvec(k);
    h01 = h0vec(k-1);
    CL1 = CLvec(k-1);
    slope = (h02 - h01)/(CL2 -CL1);
    UL = h01 + slope * (0.95 - CL1);
  end

   h0min = h0val(nh0-1);
   h0max = 5.0*h0val(nh0);

   if( CH(nh0)> 0.955)
     small = find( CH<0.945);
     if(length(small) > 0)
       h0min = h0val(length(small));
     else
       h0min = UL*0.95;
     end
       
     large = find(CH > 0.955);
     h0max = h0val(large(1));
   end

   
   fprintf(fid,' %d  %d %d \n', h0min, h0max, UL);
   bandnumber

end

fclose(fid);
