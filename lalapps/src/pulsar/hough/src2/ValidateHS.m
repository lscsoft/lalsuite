

load HS_fstatVec1.txt


Freq=HS_fstatVec1(:,1);
f1dot=HS_fstatVec1(:,4);
Alpha=HS_fstatVec1(:,2);
Delta=HS_fstatVec1(:,3);
FreqBand=0;
f1dotBand=0;
dFreq=0.00005555555555555;
df1dot=1.5432e-9;
dAlpha=0.2;
dDelta=0.2;
minStartTime=0
maxEndTime=820631477; 
refTime=600000000;
DataFiles="/local_data/badkri/fakesfts/H-1_H1*.sft";

index = 1;

cmdline = sprintf("ComputeFStatistic_v2 --Freq=%.12g --f1dot=%.12g \
--FreqBand=%.12g --f1dotBand=%.12g --dFreq=%.12g --df1dot=%.12g \
--dAlpha=%.12g  --dDelta=%.12g --minStartTime=%d --maxEndTime=%d --refTime=%d \
--DataFiles=%s --Alpha=%.12g --Delta=%.12g --AlphaBand=0.1 --DeltaBand=0.1 \
--ephemDir=/local_data/badkri/lscsoft/share/lal/ --ephemYear=05-09 \
--TwoFthreshold=0 --gridType=0 --outputLabel=CFSv2 \
--outputFstat=CFSv2 -v1", Freq(index), f1dot(index), FreqBand, f1dotBand, dFreq, df1dot, \
		  dAlpha, dDelta, minStartTime, maxEndTime, refTime, DataFiles, \
		  Alpha(index), Delta(index))  

[output, status] = system(cmdline)

load CFSv2

F_CFSv2 = CFSv2(1,7);
F_HS = 2*HS_fstatVec1(index,5);
diff=F_CFSv2 - F_HS

