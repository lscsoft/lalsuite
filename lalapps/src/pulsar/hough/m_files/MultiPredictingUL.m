%fileoutput = 'Multi_W100_1000predictedUL';
fileoutput = 'Multi_W50_100predictedUL';
fid = fopen(fileoutput, 'w');

load('MultiSigMax25_50_100W.mat');

fmin=BandList(:,1);
fmax=BandList(:,2);
signi=Sigmax';

S=signi/sqrt(2) + erfcinv(0.1);
Nbands = length(S);

kk=load('Multifit50.txt');
xx=kk(1:Nbands ,2)*2^(0.25)/1800;
x=sqrt(S).*xx;


for bandnumber = 1:Nbands;  
  fprintf(fid,'%d  %d %d %d \n', fmin(bandnumber), fmax(bandnumber),  Sigmax(bandnumber), 11.0*x(bandnumber) );
end

fclose(fid);
