
load H195UL;
load L195UL;

fMinH1 = H195UL(:,2);
fMaxH1 = H195UL(:,3);
fMinL1 = L195UL(:,2);
fMaxL1 = L195UL(:,3);


freq = 200:0.01:399.9;
nFreqs = length(freq);

for index = 1:nFreqs

  ULH1 = 1;
  ULL1 = 1;

  tempFreq = freq(index);

  temp = find(fMinH1 < tempFreq);
  if (length(temp) > 0) 
    minindex = max(temp);
    if (fMaxH1(minindex) > tempFreq)
       ULH1 = H195UL(minindex,5);
    end
  end


  temp = find(fMinL1 < tempFreq);
  if (length(temp) > 0) 
     minindex = max(temp);
    if (fMaxL1(minindex) > tempFreq)
       ULH1 = L195UL(minindex,5);
    end
  end

  UL(index) = min(ULH1, ULL1);
end

vetoindices = find(UL < 4.0*10^(-22));
freqVeto =  freq(vetoindices);
ULVeto = UL(vetoindices);
plot(freqVeto, ULVeto);
length(vetoindices)
