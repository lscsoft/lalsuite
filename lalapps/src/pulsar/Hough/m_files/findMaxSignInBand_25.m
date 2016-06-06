%kk=load('./nstar_100_1000.txt');
kk=load('./nstarMulti_50_100.txt');
significance=kk(:,2);
FreqValues=kk(:,3);

clear kk

BandList(:,1) = 50:0.25:99.75;
BandList(:,2) = BandList(:,1)+ 0.25;

Nbands = length(BandList);

for bandnumber = 1:Nbands; %the  current frequency band
   fmin = BandList(bandnumber, 1);
   fmax = BandList(bandnumber, 2);

   indices = find (FreqValues >= fmin & FreqValues <=fmax);
   if (length(indices)>0)
     Sigmax(bandnumber) = max( significance(indices));
   else
     Sigmax(bandnumber) = NaN ;
   end
   
end

save MultiSigMax25_50_100W BandList Sigmax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
