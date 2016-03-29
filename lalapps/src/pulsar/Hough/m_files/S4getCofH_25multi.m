% making up 0.25 

%prefix='/local_data/sintes/S4/MultiMC/MultiMC_100_1000_500/MCfreq_';
prefix='/local_data/sintes/S4/MultiMC/MultiMC/MCfreq_';

%fileoutput = 'MultiMC_100_1000.25_500';
fileoutput = 'MultiMC_50_100.25_500';
fid = fopen(fileoutput, 'w');

% file with driver output
%file = 'SigMax25.mat';
file = 'MultiSigMax25_50_100W.mat';
load(file);

Bands = BandList(:,1);
Nbands = length(Bands);

fshift=0;   

for bandnumber = fshift+1:(fshift+Nbands); %the  current frequency band
   fmin = BandList(bandnumber, 1);
   fmax = BandList(bandnumber, 2);
      
   basestring = strcat(prefix, num2str( fmin ) );
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
