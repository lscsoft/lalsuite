% $Id$
%
%  Matlab script to get everything
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the file with all the maximun numbercount per frequency already cleaned
%NC=load('AllSkyMax_4k_H1_0.005c');
NC=load('/home/badkri/S2results/CLEAN_L1');

%the directory where the MC results are placed
%DirectoryMC = '/local_data/badkri/S2-xxx/MC_allsky/';
%Detector = 'H1';
DirectoryMC = '/home/badkri/MCresults/';
Detector = 'L1_coarse';

fileoutput1 = strcat(Detector,'_CofH');
fid1 = fopen(fileoutput1, 'w');

fileoutput2 = strcat(Detector,'_C95');
fid2 = fopen(fileoutput2, 'w');

nh0=10; %number of different h0
Nbands = 200; %number of frequency bands

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  do not modify below this line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nc=NC(:,1);
fr=NC(:,2);
clear NC;

nc1Hz=reshape(nc,1800,Nbands);
nMax= max(nc1Hz);

freq=reshape(fr,1800,Nbands);
nFre=min(freq);

for fband=1:Nbands;
  kknc=nc1Hz(:,fband);
  Zero=find(kknc==0);
  ratio=length(Zero)/length(kknc);
  if ratio>0.9;
  	nMax(fband)=0;
  end
end


MCfilepath = strcat(DirectoryMC,Detector);
MCfilepre  = strcat(MCfilepath,'/MC_');

for bandnumber=1:Nbands;
  bn=bandnumber+199;
  
  basestring = strcat(MCfilepre, int2str( bn ) );
  ncstring = strcat(basestring, '_nc');
  h0string = strcat(basestring, '_h0');

  maxNC= nMax(bandnumber);
  UL=0.0;
  
  fprintf(fid1,'%d %d ',  nFre(bandnumber), maxNC );
  fprintf(fid2,'%d %d ',  nFre(bandnumber), maxNC );


  if maxNC>0;
     Ncount = load(ncstring);
     h0vect = load(h0string); 
     
     %getting the confidence level for each h0 value
     for h0num=1:nh0
       x=Ncount(:, h0num+1);
       kkcount = find(x>maxNC);
       CH(h0num) = length(kkcount)/length(x);
       fprintf(fid1,' %d %d ', h0vect(h0num), CH(h0num) );
     end
     fprintf(fid1,' \n');

     %%%%%%%%%%%%%%%%%%%%Refining h0
     h0min= h0vect(nh0-1);
     h0max= 5.0*h0vect(nh0);
     if( CH(nh0)> 0.96)
       small = find( CH<0.94);
       h0min = h0vect(length(small));
       large = find( CH > 0.96);
       h0max = h0vect(large(1));
     end
     fprintf(fid2,'%d %d ', h0min, h0max );

     %%%%%%%%%%%%%%%%%getting 95% upper limit
      h01=h0vect(1);
      CL1=CH(1);

      indices = find(CH > 0.95);
      if( length(indices)>0)
         k=indices(1);
         h02=h0vect(k);
         CL2=CH(k);
         if(k>1)
            h01=h0vect(k-1);
            CL1=CH(k-1);
         else
            h01 = h0vect(1);
            CL1 = CH(1);
            h02 = h0vect(2);
            CL2 = CH(2);
         end
         slope = (h02-h01)/(CL2-CL1);
         UL=h01 +slope*(0.95 -CL1);
      else
         indices = find(CH < 0.95);
         k = length(indices);
         h02 = h0vect(k);
         CL2 = CH(k);
         h01 = h0vect(k-1);
         CL1 = CH(k-1);
         slope = (h02 - h01)/(CL2 -CL1);
         UL = h01 + slope * (0.95 - CL1);
      end
      fprintf(fid2,' %d \n', UL );
 
  else
    for h0num=1:nh0
       fprintf(fid1,' 0.0 0.0 ' );   
    end
    fprintf(fid2,' 0.0 0.0 0.0 \n' );
    fprintf(fid1,' \n');
  end
 
end

 fclose(fid1);
 fclose(fid2);

