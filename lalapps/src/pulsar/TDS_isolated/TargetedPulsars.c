/*********************************************************************************/
/*                    Time domain heterodyning code for known pulsars (S2/S3)    */
/*                                                                               */
/*			               Réjean Dupuis                             */
/*                                                                               */
/*                  University of Glasgow - last modified 26/03/2004             */
/*********************************************************************************/
$Id$
#include "TargetedPulsars.h" 

#define CONDOR 1 
#define DBUG 1 /* set to 1 to run in debug mode with more print outs */

int
main( int argc, char *argv[])
{
  /***** declare variables *****/
  static LALStatus status;
  FILE *psrfp=NULL,*fpout=NULL, *fpout3=NULL;
  FILE *fpolg, *fpsf, *fpab;
  FrFile *frfile;         
  FrVect *frvect;    
  INT4 starts[MAXSTART],nstarts=0, count=0, tbase=INT_MAX, k=0;
  INT4 opencount=0, flag = 0, Nab, jobnum, i, irun, iflag;
  INT4 pnum, dec_deg, dec_min, ra_hrs, ra_min, counter = 0, i_pulsar = 0, j;
  char psrinput[256], psrout[256][MAXPSRS], chname[256];
  char framelist[256],earthfile[256],sunfile[256], outfile[256];
  char olgfile[256],sffile[256], abfile[256];
  char jnkstr1[256],jnkstr2[256],jnkstr3[256], jnkstr4[256],jnkstr5[256];   
  char pname[MAXPSRS][256], ppmra[16], ppmdec[16], pposepoch[16], ppepoch[16];  
  COMPLEX16Vector *B; 
  REAL4 wc;
  REAL8Vector *tt, *alpha, *beta, *tab;
  REAL8 RA[MAXPSRS], DEC[MAXPSRS], f0[MAXPSRS], f1[MAXPSRS], f2[MAXPSRS];
  REAL8 H0, phaseH0, C0, phaseC0, a, b, ab, tg, f=0, jnknum1;
  REAL8 ra_sec, dec_sec, pf0, pmra, pmdec,ffepoch, posepoch, fepoch_gps;
  REAL8 pf1, pf2, s_ra, s_dec;
  LIGOTimeGPS fepoch[MAXPSRS],epoch; 
  UINT4 npts;
  HeterodyneInput input;    HeterodyneOutput output;    HeterodyneParams params; 
  EphemerisData *edat = NULL;
  COMPLEX16ZPGFilter *zpgFilter = NULL;
  REAL8IIRFilter *iirFilter1Re = NULL;   REAL8IIRFilter *iirFilter1Im = NULL;
  REAL8IIRFilter *iirFilter2Re = NULL;   REAL8IIRFilter *iirFilter2Im = NULL;
  REAL8IIRFilter *iirFilter3Re = NULL;   REAL8IIRFilter *iirFilter3Im = NULL;  
  CalibrateFineParams calParams;    CalibrateFineInput calInput;  CalibrateFineOutput calOutput; 
  REAL8 REF_MJD, REF_GPS_SECS;

/* get command line arguments and print error if there is wrong number of arguments.  if 
the variable condor is defined, then read in the first argument jobnum from stdin */

#ifndef CONDOR
  if (argc !=7 || (jobnum=atoi(argv[1]))<0 || jobnum>99999){
    int a;
    fprintf(stderr,"Syntax:\n\t%s N DIR1 DIR2 DETECTOR RUN InjectionsFlag where 0<=N<=99999.\n",argv[0]);
    fprintf(stderr,"Files used are jobdata.N, jobtimes.N, and pulsar.list\n");
    fprintf(stderr,"Input files are taken from DIR1\n");
    fprintf(stderr,"Output files are put into DIR2\n");
    fprintf(stderr,"DETECTOR is either H1, H2, or L1\n");
    fprintf(stderr,"Run - 2 (S2) or 3 (S3)\n");
    fprintf(stderr,"injections? 0 for no, 1 for yes\n");
    fprintf(stderr,"There were argc=%d arguments:\n",argc);
    fflush(stderr);
    for (a=0;a<argc;a++)
	fprintf(stderr,"arg=%d: %s\n",a,argv[a]);
    fflush(stderr);
    return 2;
  }
#else
  scanf("%d",&jobnum);
  fprintf(stderr,"condor job no %d\n", jobnum);
#endif 

   irun = atoi(argv[5]);
   if (irun != 2 && irun != 3)
   {
      fprintf(stderr, "irun = ??? Should be 2 (S2) or 3 (S3)!\n");
      return(1);
   }
   
   iflag = atoi(argv[6]);
   if (iflag != 0 && iflag != 1)
   {
      fprintf(stderr, "iflag has to be 0 or 1. Either regular or hardware injections\n");
      return(0);
   }
      
/************** BEGIN READING PULSAR PARAMETERS ******************/     
/* read input file with pulsar information */
  sprintf(psrinput,"%s/pulsar.list",argv[2]);
  psrfp=tryopen(psrinput,"r");
  
  
/* if iflag = 1 then we are analysing data from a hardware injection.  the format of the input 
file is slighly different for the hardware injections */
  
  if (iflag == 1)
     fscanf(psrfp, "%d %s %lf %lf %s %s %s %lf %lf %lf %s", &pnum, &pname[counter][0], &ra_sec,&dec_sec,
     &ppmra[0], &ppmdec[0], &pposepoch[0], &pf0, &pf1, &pf2,  &ppepoch[0]); 
  else 
     fscanf(psrfp, "%d%s %d:%d:%lf %d:%d:%lf %s %s %s %lf %lf %lf %s",
     &pnum, &pname[counter][0], &ra_hrs, &ra_min, &ra_sec, &dec_deg, &dec_min, &dec_sec,
     &ppmra[0], &ppmdec[0], &pposepoch[0], &pf0, &pf1, &pf2,  &ppepoch[0]); 
     
  while (!feof(psrfp))
  {
    f0[counter] = pf0;
    f1[counter] = pf1;
    f2[counter] = pf2;
    
    /* if this is a hardware injection then divide the spin parameters provided by 2.0 since
    the rest of the code is especting the pulsar rotation frequency as opposed to the gw frequency */
    if (iflag == 1){
      f0[counter] /= 2.0;
      f1[counter] /= 2.0;
      f2[counter] /= 2.0;
    }
    if ((2.0*f0[counter]) > 1500.0 || (2.0*f0[counter]) < 10.0)
    {
       fprintf(stderr, "f0 outside range for truncated calibration info!\n");
       return(6);
    }
   
    /* determine sign of RA and DEC 
      this will fail if ra_hrs or dec_deg 
      are -0 (meaning that min, sec are neg)
      there aren't currently any pulsar that fit that criteria, 
      but will need to improve this parsing when we do */
        
    if (ra_hrs < 0) s_ra = -1.0;
    else s_ra = 1.0;
    if (ra_hrs == 0) fprintf(stderr, "Warning: ra_hrs = 0 - make sure conversion to radians correct\n");
    
    if (dec_deg < 0) s_dec = -1.0;
    else s_dec = 1.0;
    if (dec_deg == 0) fprintf(stderr, "Warning: dec_deg = 0 - make sure conversion to radians correct\n");

    /* again the input files for hardware injections provide info in different format
    than regular pulsars (RA and DEC in radians directly as opposed to hours (or degress) minutes, secs */
    if (iflag ==1 ){
      RA[counter] = ra_sec;
      DEC[counter] = dec_sec;}
    else {
      RA[counter] = 2.0*LAL_PI*((ra_hrs + s_ra*ra_min/60.0 + s_ra*ra_sec/3600.0)/24.0);
      DEC[counter] = 2.0*LAL_PI*((dec_deg + s_dec*dec_min/60.0 + s_dec*dec_sec/3600.0)/360.0);}
    
    if (strcmp(ppmra, "*")) { pmra = atof(ppmra);}
      else {pmra = 0.0;}
    if (strcmp(ppmdec, "*")) { pmdec = atof(ppmdec);}
      else {pmdec = 0.0;}
    if (strcmp(pposepoch, "*")) { posepoch = atof(pposepoch);}
      else {posepoch = 0.0;}
    if (strcmp(ppepoch, "*")) { ffepoch = atof(ppepoch);}
      else {ffepoch = 0.0;}
  
    /*adjust RA and DEC to epoch at beg of run */
    if (irun==2)
    {
       REF_MJD = S2_MJD;
       REF_GPS_SECS = S2_GPS_SECS;
    }
    else if (irun==3)
    {
       REF_MJD = S3_MJD;
       REF_GPS_SECS = S3_GPS_SECS;
    }
    else 
    {
      fprintf(stderr, "irun must be either 2 or 3!\n");
      return(1);
    }   
    
    if (iflag != 1) { 
      RA[counter] +=  ((REF_MJD - posepoch)/365.26) * (2.0*LAL_PI*pmra/(3600.0*360.0*1000.0));
      DEC[counter] += ((REF_MJD - posepoch)/365.26) * (2.0*LAL_PI*pmdec/(3600.0*360.0*1000.0));}
 
    posepoch = REF_MJD;
    if (ffepoch != 0.0) fepoch_gps = (-REF_MJD + ffepoch)*86400.0 + REF_GPS_SECS; 
    else fepoch_gps = REF_GPS_SECS;
   
    fepoch[counter].gpsSeconds = (INT4)floor(fepoch_gps);
    fepoch[counter].gpsNanoSeconds = (INT4)floor((fmod(fepoch_gps, 1.0)*1.e9)); 

/* if S2 hardware injections then hardwire the frequency epochs here
- due to fact that I did not have these epochs in MJD format so easier to 
just put them in manually here */	    
    if (iflag==1 && irun==2)
    {      
      if ((f0[counter]*2.0) < 1285.0) /* thus signal 1*/
      {
        fepoch[counter].gpsSeconds =733967667;
        fepoch[counter].gpsNanoSeconds = 26112310;
      }
     else if ((f0[counter]*2.0) > 1285.0) /* thus signal 2*/
     {
        fepoch[counter].gpsSeconds =733967751;
        fepoch[counter].gpsNanoSeconds =  522490380;
     }
     else
     {fprintf(stderr, "frequency of signal should not be exactly 1285Hz!\n"); return(0);}    
    }
  
    if (DBUG) fprintf(stderr, "%s\tfepoch.gpsSeconds=%d\nfepoch.gpsNanoSeconds=%d\n", 
        pname[counter], fepoch[counter].gpsSeconds, fepoch[counter].gpsNanoSeconds);

    
    /* print out pulsar information read to output directory, this file will be used later
    in the ComputePDF code */
    sprintf(psrout[counter],"%s/%s/%s/%s",argv[3],argv[4], pname[counter],pname[counter]);
    fpout3=tryopen(psrout[counter],"w");
    fprintf(fpout3, "ra\t%f\ndec\t%f\npmra\t%f\npmdec\t%f\npepoch\t%f\nf0\t%f\nf1\t%e\nf2\t%e\nfepoch\t%f",
		     RA[counter], DEC[counter], pmra, pmdec, posepoch, f0[counter],
		     f1[counter], f2[counter], ffepoch); 
    fflush(fpout3);
    fclose(fpout3);
    counter++;
   
   /* read another line of the pulsar input file */
    if (!(feof(psrfp))) 
    {
      if (iflag==1) 
          fscanf(psrfp, "%d %s %lf %lf %s %s %s %lf %lf %lf %s",  &pnum, &pname[counter][0], &ra_sec,&dec_sec,
          &ppmra[0], &ppmdec[0], &pposepoch[0], &pf0, &pf1, &pf2,  &ppepoch[0]); 
      else
          fscanf(psrfp, "%d%s %d:%d:%lf %d:%d:%lf %s %s %s %lf %lf %lf %s",
          &pnum, &pname[counter][0], &ra_hrs, &ra_min, &ra_sec, &dec_deg, &dec_min, &dec_sec,
          &ppmra[0], &ppmdec[0], &pposepoch[0], &pf0, &pf1, &pf2,  &ppepoch[0]); 
    }
  }  /* loop until all the lines are read from the pulsar.list file */
  fclose(psrfp);
  /************** END READING PULSAR PARAMETERS ******************/   
  
  /************** INITIALIZE 'STATIC' VARIABLES BEFORE GOING INTO LOOP OVER PULSARS *******/
  if (!strcmp(argv[4],"H1") || !strcmp(argv[4],"H2")) 
    params.detector = lalCachedDetectors[LALDetectorIndexLHODIFF];      
  else if (   !strcmp(argv[4],"L1"))
   params.detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  else if (   !strcmp(argv[4],"GEO"))
   params.detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];             
  else
  {
    fprintf(stderr,"Error DETECTOR must be either H1, H1, L1, or GEO\n");
    return 2;
  }  
    
  /* construct frame channel name */   
  if (!strcmp(argv[4],"GEO"))
   sprintf(chname,"G1:DER_DATA_H");
  else 
    sprintf(chname,"%s:LSC-AS_Q",argv[4]);
   
   npts = SRATE*60;
   
   /* allocate memory for 60 second chucks of raw data*/
   input.V.data = NULL;
   LALDCreateVector( &status, &input.V.data, npts);  
   input.V.deltaT = 1.0/SRATE;
   input.V.data->length = npts; 

   /* allocate memory for set of Bk's */
   B = NULL;
   LALZCreateVector( &status, &B, MAXFINE);
   TESTSTATUS(&status);
  
  /* allocate memory for time stamps for each Bk */
   tt = NULL;
   LALDCreateVector( &status, &tt, MAXFINE);          
   TESTSTATUS(&status);


  /* allocate memory to store the alpha/beta calibration parameters */   
   alpha = NULL;
   LALDCreateVector(&status, &alpha, MAXLINES);
   TESTSTATUS(&status);
  
   beta = NULL;
   LALDCreateVector(&status, &beta, MAXLINES);
   TESTSTATUS(&status);
   
  /* allocate memory for time stamps for each alpha/beta */    
   tab = NULL;
   LALDCreateVector(&status, &tab, MAXLINES);
   TESTSTATUS(&status);  
      
  /* set up three Butterworth 3 order filters before looping over pulsars 
      cut off frequencies fc1, fc2, and fc3 are defined in TargetedPulsars.h */
   wc = tan(LAL_PI * fc1 / SRATE);
   LALCreateCOMPLEX16ZPGFilter( &status, &zpgFilter, 0, 3 );
   zpgFilter->poles->data[0].re = wc*SQRT3_2;
   zpgFilter->poles->data[0].im = wc*0.5;
   zpgFilter->poles->data[1].re = 0.0;
   zpgFilter->poles->data[1].im = wc;
   zpgFilter->poles->data[2].re = -wc*SQRT3_2;
   zpgFilter->poles->data[2].im = wc*0.5;
   zpgFilter->gain.re = 0.0;
   zpgFilter->gain.im = wc*wc*wc;
   LALWToZCOMPLEX16ZPGFilter( &status, zpgFilter );
   LALCreateREAL8IIRFilter( &status, &iirFilter1Re, zpgFilter );
   LALCreateREAL8IIRFilter( &status, &iirFilter1Im, zpgFilter );
   LALDestroyCOMPLEX16ZPGFilter( &status, &zpgFilter );
   params.iirFilter1Re = iirFilter1Re;
   params.iirFilter1Im = iirFilter1Im;
    
   wc = tan(LAL_PI * fc2 / SRATE);
   LALCreateCOMPLEX16ZPGFilter( &status, &zpgFilter, 0, 3 );
   zpgFilter->poles->data[0].re = wc*SQRT3_2;
   zpgFilter->poles->data[0].im = wc*0.5;
   zpgFilter->poles->data[1].re = 0.0;
   zpgFilter->poles->data[1].im = wc;
   zpgFilter->poles->data[2].re = -wc*SQRT3_2;
   zpgFilter->poles->data[2].im = wc*0.5;
   zpgFilter->gain.re = 0.0;
   zpgFilter->gain.im = wc*wc*wc;
   LALWToZCOMPLEX16ZPGFilter( &status, zpgFilter );
   LALCreateREAL8IIRFilter( &status, &iirFilter2Re, zpgFilter );
   LALCreateREAL8IIRFilter( &status, &iirFilter2Im, zpgFilter );
   LALDestroyCOMPLEX16ZPGFilter( &status, &zpgFilter );
   params.iirFilter2Re = iirFilter2Re;
   params.iirFilter2Im = iirFilter2Im;  
   
   wc = tan(LAL_PI * fc3 / SRATE);
   LALCreateCOMPLEX16ZPGFilter( &status, &zpgFilter, 0, 3 );
   zpgFilter->poles->data[0].re = wc*SQRT3_2;
   zpgFilter->poles->data[0].im = wc*0.5;
   zpgFilter->poles->data[1].re = 0.0;
   zpgFilter->poles->data[1].im = wc;
   zpgFilter->poles->data[2].re = -wc*SQRT3_2;
   zpgFilter->poles->data[2].im = wc*0.5;
   zpgFilter->gain.re = 0.0;
   zpgFilter->gain.im = wc*wc*wc;
   LALWToZCOMPLEX16ZPGFilter( &status, zpgFilter );
   LALCreateREAL8IIRFilter( &status, &iirFilter3Re, zpgFilter );
   LALCreateREAL8IIRFilter( &status, &iirFilter3Im, zpgFilter );
   LALDestroyCOMPLEX16ZPGFilter( &status, &zpgFilter );
   params.iirFilter3Re = iirFilter3Re;
   params.iirFilter3Im = iirFilter3Im;  
  
  /*********** begin loop over pulsars ***************/ 
  for (i_pulsar=0;i_pulsar<counter;i_pulsar++)  
  { 
 /* counts number of minutes processed (or number of Bk's produced) */
    k = 0; 
  
 /* set up input paramters for LALHeterodyneToPulsar for pulsar i_pulsar */ 
   input.V.f0 = 0.0; /* base frequency */
   input.f0 = 2.0*f0[i_pulsar];
   input.f1 = 2.0*f1[i_pulsar];
   input.f2 = 2.0*f2[i_pulsar]; 
   input.source.longitude = RA[i_pulsar];
   input.source.latitude = DEC[i_pulsar];
   input.fEpochGPS = (double)fepoch[i_pulsar].gpsSeconds + (double)fepoch[i_pulsar].gpsNanoSeconds*1.0e-9;
   input.V.deltaT = 1.0/SRATE; 
   
   if (DBUG) fprintf(stderr, "%s\tf0 = %f\tf1 = %e\tf2= %e\n", pname[i_pulsar], f0[i_pulsar], f1[i_pulsar], f2[i_pulsar]);  
 
 /* begin reading calibration information for this pulsar */   
   if (strcmp(argv[4],"GEO"))
    {
      sprintf(olgfile,"%s/calibration/%sopenloopgain.txt",argv[2],argv[4]);    
      fpolg = tryopen(olgfile,"r");
    
      f = 0.0;
      while(f<2.0*f0[i_pulsar])  
         fscanf(fpolg,"%lf%lf%lf",&f,&H0,&phaseH0);       
       
      if (fabs(f-2.0*f0[i_pulsar]) > 0.1) 
         {fprintf(stderr,"Open loop gain file does not have frequency bin near f0\n");  return(4);}
      fclose(fpolg);    
      
      sprintf(sffile,"%s/calibration/%sASQsensing.txt",argv[2],argv[4]);   
      fpsf = tryopen(sffile,"r");    
      f = 0.0;
      while(f<2.0*f0[i_pulsar])
         fscanf(fpsf,"%lf%lf%lf",&f,&C0,&phaseC0);
    
      if (fabs(f-2.0*f0[i_pulsar]) > 0.1) 
         {fprintf(stderr,"Sensing function file does not have frequency bin near f0\n"); return(4);}
      fclose(fpsf);         
    
      calParams.H0 = H0; 
      calParams.phaseH0 = phaseH0;
      calParams.C0 = C0; 
      calParams.phaseC0 = phaseC0;
      
      sprintf(abfile,"%s/calibration/%salphabeta.txt",argv[2],argv[4]);
  
      fpab = tryopen(abfile,"r");
      fscanf(fpab, "%s%s%s%s%s", &jnkstr1[0], &jnkstr2[0], &jnkstr3[0], &jnkstr4[0], &jnkstr5[0]);
    
      i = 0;
      while (!feof(fpab)){
        fscanf(fpab, "%lf%lf%lf%lf%lf", &tg, &ab, &a, &b, &jnknum1);
        if (a>=ALPHAMIN && a<=ALPHAMAX){
          tab->data[i] = tg;
          alpha->data[i] = a;
          beta->data[i] = b;
          i++;
  	  if (i > MAXLINES){ fprintf(stderr,"MAXLINES not big enough!\n");return(4);}
        }
      }
      fclose(fpab);
      Nab = i;
    } 
  /* end reading calibration information for this pulsar */  
    
    /* beg read frame data (only if first pulsar) */
    if (i_pulsar==0)
    {   
      {
       char timelist[256];
       INT4 tbasei;
       FILE *stfp=NULL;
       sprintf(timelist,"%s/data%s/jobtimes.%05d",argv[2],argv[4],jobnum);
       stfp=tryopen(timelist,"r"); 
       while (2==fscanf(stfp,"%d %d",starts+nstarts,&tbasei)){
         if (nstarts==0)
           tbase=tbasei;
         else if (tbasei!=tbase){
           fprintf(stderr,"File %s contains inconsistent SFT time baselines %d != %d\n",timelist,tbasei,tbase);
	   fflush(stderr);
	   return 3;
         }
         nstarts++;
         if (nstarts>=MAXSTART){
	   fprintf(stderr,"More than MAXSTART=%d lines in file: %s.  Increase MAXSTART\n",
	   MAXSTART,timelist);
	   fflush(stderr);
	   return 3;
         }
       }
       if (!nstarts){
         fprintf(stderr,"File %s didn't contain any valid lines!\n",timelist);
         fflush(stderr);
         return 3;
       }
       fclose(stfp);
      }  
      tbase = TBASE;
      npts = SRATE*tbase; 
      sprintf(framelist,"%s/data%s/jobdata.%05d.ffl",argv[2],argv[4],jobnum);
      opencount=0;
      while (!(frfile = FrFileINew(framelist))){
        fprintf(stderr, "Couldnt open frame file %s\n", framelist);
        fflush(stderr);
      if (opencount++<10)
        sleep(10);
      else
        return(3);
      }
    } 
     /* end read frame data (for first pulsar) */

 /* beg extracting from frvect and heterodyning the data*/
  for (count=0;count<nstarts;count++) 
  {
    epoch.gpsSeconds=starts[count];
    epoch.gpsNanoSeconds=0;    
    input.V.epoch = epoch;
    if (DBUG) 
      fprintf(stderr, "%d. %s -- %d out of %d\n", i_pulsar,pname[i_pulsar], count, nstarts);fflush(stderr);
 
    frvect = FrFileIGetVAdc(frfile, chname, epoch.gpsSeconds, tbase, 0);
    if (frvect == NULL)
    {
      fprintf(stderr, "Data between times %d and %d not present\n",epoch.gpsSeconds,epoch.gpsSeconds+tbase);
      fflush(stderr);
      continue;
    }
    if (strcmp(argv[4],"GEO"))
    {     
      if ( frvect->type != FR_VECT_4R )
      {
        fprintf(stderr, "Wrong data type (not FR_VECT_4R) found in frame!\n" );
        fflush(stderr);
        return(5);
      }   
    } 
    else /* if GEO */
    {
      if ( frvect->type != FR_VECT_8R )
      {
        fprintf(stderr, "Wrong data type (not FR_VECT_8R) found in frame!\n" );
        fflush(stderr);
        return(5);
      }   
    }
    
    if (frvect->next)
   {
       fprintf(stderr, "Data between times %d and %d had a gap\n",epoch.gpsSeconds,epoch.gpsSeconds+tbase);
       fflush(stderr);
       FrVectFree(frvect);
       frvect=NULL;
       continue;
   }    
   
   if (strcmp(argv[4],"GEO"))  /* if LIGO */
   {     
     for (i=0;i<npts;i++) 
       input.V.data->data[i] = frvect->dataF[i];
   }
   else /* if GEO */
   {
     for (i=0;i<npts;i++) 
       input.V.data->data[i] = frvect->dataD[i];
   }
   
   FrVectFree(frvect);
   frvect=NULL;

  /* set up LALBarycenter */
   edat = (EphemerisData *)LALMalloc(sizeof(EphemerisData));    
   (*edat).leap = 13; 
   sprintf(earthfile,"%s/earth00-04.dat",argv[2]);
   sprintf(sunfile,"%s/sun00-04.dat",argv[2]);
   (*edat).ephiles.earthEphemeris = earthfile;
   (*edat).ephiles.sunEphemeris = sunfile;	 
   LALInitBarycenter(&status, edat);
   if(status.statusCode) 
   {
     fprintf(stderr, "Unexpectedly got error code %d and message %s\n", 
     status.statusCode, status.statusDescription);
     return 0;
   } 
   params.edat = edat;

/******* CALL LALHeterodyneToPulsar ************/     
   LALHeterodyneToPulsar( &status, &output, &input, &params ); 
   if(status.statusCode) 
   {
     fprintf(stderr, "Unexpectedly got error code %d and message %s\n", 
     status.statusCode, status.statusDescription);
     return 0;
   }      

/* store Bk's and time stamps into vectors */
   if (k>0)/* due to  to iir step response do not write first minute after start */
   {   
       B->data[k-1].re =  output.B.re;
       B->data[k-1].im =  output.B.im; 
       tt->data[k-1] = (double)epoch.gpsSeconds + (double)epoch.gpsNanoSeconds*1.0e-9;  
   }

   k++; 
   LALFree(edat->ephemE);
   LALFree(edat->ephemS);
   LALFree(edat);
  }   /***** end extracting from frvect and heterodyning the data *****/

 /**** begin calibrating Bk's and writing Bk's to file****/   
 if (DBUG) fprintf(stderr, "opening output file: epoch %d, psrname %s\n", epoch.gpsSeconds, pname[i_pulsar]); 
 
 /* create output file name and try to open (will fail if output director doesn't exist) */
 sprintf(outfile,"%s/%s/%s/outfine.%s_%s.%d",argv[3],argv[4],pname[i_pulsar],pname[i_pulsar],argv[4],epoch.gpsSeconds);
 fpout=tryopen(outfile,"w");    
 
 flag=0;
 /* only calibrate if we have LIGO data */
 if (strcmp(argv[4],"GEO"))
  {   
    /* for each Bk (total of k-1 Bk's) first calibrate then write out */
    for (i=0;i<k-1;i++) { 
      /* calibration information has already been read. however, there are two reference 
      	 epochs for H2 S2.  so if this is H2 S2 and we need to use the second epoch, then
	 read the second reference ASQ sensing and open loop gain */ 
      if (!strcmp(argv[4],"H2") && tt->data[i] > 731849043.0 && flag == 0 && irun == 2)
      {
        flag = 1;
	if (DBUG) fprintf(stderr, "reading second epoch of S2 calibration for H2\n");
        sprintf(olgfile,"%s/calibration/%sopenloopgainb.txt",argv[2],argv[4]);
        fpolg = tryopen(olgfile,"r");
        f = 0.0;
        while(f<2.0*f0[i_pulsar])  
         fscanf(fpolg,"%lf%lf%lf",&f,&H0,&phaseH0);
        if (fabs(f-2.0*f0[i_pulsar]) > 0.1) 
         {fprintf(stderr,"Open loop gain file does not have frequency bin near f0\n");  return(4);}
        fclose(fpolg);    
        sprintf(sffile,"%s/calibration/%sASQsensingb.txt",argv[2],argv[4]);  
        fpsf = tryopen(sffile,"r");    
        f = 0.0;
        while(f<2.0*f0[i_pulsar])
         fscanf(fpsf,"%lf%lf%lf",&f,&C0,&phaseC0);
        if (fabs(f-2.0*f0[i_pulsar]) > 0.1) 
       {fprintf(stderr,"Sensing function file does not have frequency bin near f0\n");return(4);}
       fclose(fpsf);         
       calParams.H0 = H0; 
       calParams.phaseH0 = phaseH0;
       calParams.C0 = C0; 
       calParams.phaseC0 = phaseC0;       
    } 
      /* find first alpha/beta time stamp that is less than 60 seconds away from Bk time stamp */    
      for (j=0;j<Nab;j++) {
        if(fabs(tt->data[i] - tab->data[j])<=30)  
        {
  	  calInput.alpha = alpha->data[j];
	  calInput.beta = beta->data[j];
	  calInput.B.re = B->data[i].re;
	  calInput.B.im = B->data[i].im;
	  
          /* calibrate Bk's and write to file */
	  LALCalibrateFineHeterodyne(&status, &calOutput, &calInput, &calParams);
	  fprintf(fpout,"%f\t%e\t%e\n", tt->data[i], calOutput.B.re, calOutput.B.im);
          break;
        }
      } /* loop over alpha/beta */
    } /* loop over all Bk's */
  } /* if LIGO data then calibrate */
  
  else /* if GEO no need to calibrate just write to file */
  {
    for (i=0;i<k-1;i++)  
       fprintf(fpout,"%f\t%e\t%e\n", tt->data[i], B->data[i].re, B->data[i].im);              
  }  
  
  fclose(fpout);
 /**** end calibrating Bk's and writing Bk's to file****/   
  flag = 0;
} /*********** end loop over pulsars ***************/ 
  
  /* free memory and check for memory leaks */
  LALZDestroyVector(&status, &B);   
  LALDDestroyVector(&status, &tt);
  LALDDestroyVector(&status, &tab);   
  LALDDestroyVector(&status, &alpha);   
  LALDDestroyVector(&status, &beta);   
  LALDestroyREAL8IIRFilter( &status, &iirFilter1Re );
  LALDestroyREAL8IIRFilter( &status, &iirFilter1Im ); 
  LALDestroyREAL8IIRFilter( &status, &iirFilter2Re );
  LALDestroyREAL8IIRFilter( &status, &iirFilter2Im );   
  LALDestroyREAL8IIRFilter( &status, &iirFilter3Re );
  LALDestroyREAL8IIRFilter( &status, &iirFilter3Im );  
  LALDDestroyVector(&status, &input.V.data); 
  TESTSTATUS( &status );
  LALCheckMemoryLeaks();
  return 0;
} /* main */
