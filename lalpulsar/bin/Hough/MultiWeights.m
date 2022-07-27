#! /usr/bin/octave

## Dumps weights vectors or noise estimate for a set of sfts

## ./MultiWeights usage
## 
##  -d                    INT      set lalDebugLevel [0] 
##   -h, --help            BOOL     Print this message [] 
##   -f, --f0              REAL     Start search frequency [505] 
##   -b, --fSearchBand     REAL     Search frequency band [0.05] 
##       --printLog        BOOL     Print Log file [FALSE] 
##   -D, --sftDir          STRING   SFT filename pattern [REQUIRED] 
##       --linefiles       LIST     Comma separated List of linefiles (filenames must contain IFO name) [NULL] 
##       --startTime       REAL     GPS start time of observation [0.0] 
##       --endTime         REAL     GPS end time of observation [2147483647] 
##       --timeStampsFile  STRING   Input time-stamps file [NULL] 
##       --weightAM        BOOL     Use amplitude modulation weights [TRUE] 
##       --weightNoise     BOOL     Use SFT noise weights [TRUE] 
##       --dumpAllWeights  BOOL     Dump all weights [FALSE] 
##       --dumpRelativeWeights BOOL     Dump IFO relative weights [FALSE] 
##       --dumpNoise       BOOL     Dump Noise estimate [FALSE] 
##   -E, --earthEphemeris  STRING   Earth Ephemeris file ["./earth05-09.dat"] 
##   -S, --sunEphemeris    STRING   Sun Ephemeris file ["./sun05-09.dat"] 
##       --AlphaWeight     REAL     sky Alpha for weight calculation [1] 
##       --DeltaWeight     REAL     sky Delta for weight calculation [1] 
##       --outfile         STRING   output file name ["./tempout"] 


## DataFiles="/sft/S5-LIGO/H1-1800SFT/*.sft";
## DataFiles="/sft/S5-LIGO/H2-1800SFT/*.sft";
## DataFiles="/sft/S5-LIGO/L1-1800SFT/*.sft";

DataFiles="/sft/S5-links/*.sft";
EarthE="/home/sintes/earth05-09.dat";
SunE="/home/sintes/sun05-09.dat";

## startT=815410983.0;
## stopT=816335150.0;

DumpAllWeights=1;
DumpRelativeW=1;
DumpNoise=0;

WeightAM=1;
WeightNoise=1;

fStart=93.5;
fBand=0.05;

AlphaWeight=0.0;
DeltaWeight=1.0;

system("mkdir -p ./MultiWeightsDump");

cmdline = sprintf("/home/all64/lscsoft/lalapps/src/pulsar/hough/src/MultiWeights  --sftDir='%s' \
         --earthEphemeris=%s  --sunEphemeris=%s \
         --f0=%.12g --fSearchBand=%.12g  \
	 --AlphaWeight=%.12g --DeltaWeight=%.12g \
	 --weightAM=%.12g --weightNoise=%.12g \
	 --dumpAllWeights=%.12g --dumpRelativeWeights=%.12g --dumpNoise=%.12g -d 0",\
       DataFiles, EarthE, SunE, fStart, fBand,  AlphaWeight, DeltaWeight,\
       WeightAM, WeightNoise, DumpAllWeights, DumpRelativeW,DumpNoise); 

system ( cmdline );

cat = "cat ./tempout >> ./MultiWeightsDump/MultiWeights_f93.5_N1AM1_delta1.0.txt";
system (cat);
system("rm -f ./tempout");



WeightAM=0;

cmdline = sprintf("/home/all64/lscsoft/lalapps/src/pulsar/hough/src/MultiWeights  --sftDir='%s' \
         --earthEphemeris=%s  --sunEphemeris=%s \
         --f0=%.12g --fSearchBand=%.12g  \
	 --AlphaWeight=%.12g --DeltaWeight=%.12g \
	 --weightAM=%.12g --weightNoise=%.12g \
	 --dumpAllWeights=%.12g --dumpRelativeWeights=%.12g --dumpNoise=%.12g -d 0",\
       DataFiles, EarthE, SunE, fStart, fBand,  AlphaWeight, DeltaWeight,\
       WeightAM, WeightNoise, DumpAllWeights, DumpRelativeW,DumpNoise); 

system ( cmdline );

cat = "cat ./tempout >> ./MultiWeightsDump/MultiWeights_f93.5_N1AM0.txt";
system (cat);
system("rm -f ./tempout");
