kk=load('./nstar_800_900_noweights.txt');
significance=kk(:,2);
FreqValues=kk(:,3);

clear kk

BandList(:,1) = 800:0.5:899.5;
BandList(:,2) = BandList(:,1)+ 0.5;

Nbands = length(BandList);

for bandnumber = 1:Nbands; %the  current frequency band
   fmin = BandList(bandnumber, 1);
   fmax = BandList(bandnumber, 2);

   indices = find (FreqValues >= fmin & FreqValues <=fmax);
   Sigmax(bandnumber) = max( significance(indices));
end

save SigMaxNoW800 BandList Sigmax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
