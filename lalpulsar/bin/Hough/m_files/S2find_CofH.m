%
%  Matlab script to get everything
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the file with all the maximun numbercount per frequency already cleaned
%NC=load('AllSkyMax_4k_H1_0.005c');
%NC=load('/local_data/badkri/S2results/CLEAN_H1');

NC=load('/home/sintes/S2_ligo/AllSkyMax_2k_H2_0.005c2');
%NC=load('/home/sintes/S2_ligo/AllSkyMax_4k_H1_0.005c2');
%NC=load('/home/sintes/S2_ligo/AllSkyMax_4k_L1_0.005c');

%the directory where the MC results are placed
DirectoryMC = '/home/sintes/S2_ligo/MC_allsky/';
Detector = 'H2';
%Detector = 'H1_newline_finest';
%Detector = 'L1_fine';

MCfilepath = strcat(DirectoryMC,Detector);
MCfilepre  = strcat(MCfilepath,'/MC_H2_');
%MCfilepre  = strcat(MCfilepath,'/MC_H1_');
%MCfilepre  = strcat(MCfilepath,'/MC_L1_');

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
       kkcount = find(x>=maxNC);
       CH(h0num) = length(kkcount)/length(x);
       fprintf(fid1,' %d %d ', h0vect(h0num), CH(h0num) );
     end
     fprintf(fid1,' \n');

     %%%%%%%%%%%%%%%%%%%%Refining h0
     h0min= h0vect(nh0-1);
     h0max= h0vect(nh0)/CH(nh0);
     if( CH(nh0)> 0.955)
       small = find( CH<0.945);
       if (length(small) ~= 0)
              h0min = h0vect(length(small));
       elseif (CH(1) < 0.95)
              h0min = h0vect(1);
       else
              h0min = 0.95*h0vect(1);
       end
       large = find( CH > 0.955);
       h0max = h0vect(large(1));
     end


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
         ii = 1;
         indices = find(CH < 0.95);
         k = length(indices);
         h02 = h0vect(k);
         CL2 = CH(k);
         h01 = h0vect(k-ii);
         CL1 = CH(k-ii);
         while (CL1 == CL2)
              ii = ii+1;
              CL1 = CH(k-ii);
              h01 = h0vect(k-ii);
         end
         slope = (h02 - h01)/(CL2 -CL1);
         UL = h01 + slope * (0.95 - CL1);
      end

      if (UL < h0min)
         h0min = 0.97 * UL;
      end
      if (UL > h0max)
         h0max = 1.03 * UL;
      end
      fprintf(fid2,'%d %d %d\n', h0min, h0max, UL );

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
