kk=load('./nstar_100_1000.txt');
significance=kk(:,2);
FreqValues=kk(:,3);

clear kk

BandList(:,1) = 100:0.5:999.5;
BandList(:,2) = BandList(:,1)+ 0.5;

Nbands = length(BandList);

for bandnumber = 1:Nbands; %the  current frequency band
   fmin = BandList(bandnumber, 1);
   fmax = BandList(bandnumber, 2);

   indices = find (FreqValues >= fmin & FreqValues <=fmax);
   Sigmax(bandnumber) = max( significance(indices));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
