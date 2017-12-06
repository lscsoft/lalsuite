/*
*  Copyright (C) 2007 Reinhard Prix, Yousuke Itoh
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/**
 * \author Yousuke Itoh
 */

/* 
  gcc -Wall -g -o FstatShapeTest FstatShapeTest.c -lm -I/afs/aeiw/grawave/Linux/lal/include/ -llal -llalsupport -L/afs/aeiw/grawave/Linux/lal/lib/
*/
/*
  ./FstatShapeTest -o FaFb00.001 -t FaFb01.001 > FST.txt
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>

#include "getopt.h"

#define MAX_DATA_POINTS 100000000

/* #define DEBG_CONSTDATA  */

/* function prototype */ 
INT2 InitData(void);
INT2 HandleComArg(INT4 argc, CHAR ** argv);
INT2 ReadHeader(void);
INT2 ReadData(void);
INT2 SummitFinder(void);
INT2 ShiftData(REAL8);
INT2 RearrangeData(void);
INT2 ConstructData(void);
INT2 ComputeChi(void);
REAL8 chi2cdf(REAL8, REAL8);
INT4 myRound(REAL8);
REAL8 gammaq(REAL8 dof, REAL8 chi2var);
void gcf(REAL8 * gammacf, REAL8 a, REAL8 x, REAL8 * gln);
void gser(REAL8 * gamser, REAL8 a, REAL8 x, REAL8 * gln);
REAL8 gammln(REAL8 xx);

/* structure typedef */
typedef struct tagGlobalVariables {
  INT4  nData;
  INT4  nSearch;
  INT4  indexOffset;
  INT4  indexStepO;
  INT4  indexStepT;
  REAL8 sFreq;   /* starting frequency of the registered cluster */
  REAL8 dFreq;   /* frequency resolution */
  REAL8 fmax;    /* the frequency where the maximum amplitude occurs */
  REAL8 Acoeff;  /* A coefficient of JKS = (a(t)||a(t)) */
  REAL8 Bcoeff;  /* B coefficient of JKS = (b(t)||b(t)) */
  REAL8 Ccoeff;  /* C coefficient of JKS = (a(t)||b(t)) */
  REAL8 * lmfrq; 
} GlobalVariables;

typedef struct tagHeader {
  INT4  nData;
  REAL8 sFreq;
  REAL8 dFreq;
  REAL8 fmax;
  REAL8 FSmax;
} Header;

typedef struct tagFaFb {
  REAL8 freq;
  REAL8 RFa;
  REAL8 IFa;
  REAL8 RFb;
  REAL8 IFb;
  REAL8 Fstat;
} FaFbdata;


/* variables declaration */
GlobalVariables GV;
Header ObsvHeader, TestHeader;
FaFbdata *FaFbObsv, *FaFbTest,*FaFbSubt;

REAL8 ChistatSum=0.0;     /* Chi square statistic  */
REAL8 sigLevel=0.01;

/* static LALStatus status; */

/* default input file name */
const CHAR * obsvdatafile = "FaFbObsv.txt";
const CHAR * testdatafile = "FaFbTest.txt";

/* flag */
INT2 verboseflag=0; /* verbose flag */
INT2 computeProbFlag=0; 

/*-----------------------------------------------------------------*/
/*                                                                 */
/*-----------------------------------------------------------------*/


INT4 main(INT4 argc, CHAR ** argv)
{
  INT4 iiter;
  REAL8 lmfrq,fmaxT,sFreqT;

  REAL8 df=0.0;     /* degrees of freedom of the chi square distribution */

  REAL8 ChiSumMin=FLT_MAX;
  REAL8 lmfrqmin=0.0;

  REAL8 prob;
  INT2 rejection=0;

  if(InitData()) return 1;

  if(HandleComArg(argc, argv)) return 11;

  if(ReadHeader()) return 31;

  if(ReadData()) return 41;

  if(SummitFinder()) return 51;


  fmaxT =TestHeader.fmax;
  sFreqT=TestHeader.sFreq;

  /* shift the data */
  for(iiter=0;iiter<GV.nSearch;iiter++) {

    TestHeader.fmax=fmaxT;
    TestHeader.sFreq=sFreqT;
    lmfrq=GV.lmfrq[iiter];

    if(ShiftData(lmfrq))
      return 61;

  /* rearrange the data if necessary */
    if(
       ((ObsvHeader.nData)!=(TestHeader.nData))||
       ((ObsvHeader.sFreq)!=(TestHeader.sFreq))||
       ((ObsvHeader.dFreq)!=(TestHeader.dFreq))
       ) {
      if(RearrangeData())  
	return 71;
    } else {
      GV.nData=ObsvHeader.nData;
      GV.sFreq=ObsvHeader.sFreq;
      GV.dFreq=ObsvHeader.dFreq;
      GV.indexOffset=0;
      GV.indexStepO=1;
      GV.indexStepT=1;
    }

    if(ConstructData()) return 81;


    if(ComputeChi()) return 91;


  /* chi square test */
  /* \sum_{n=0}^{nData}2 F[irec] is our statistic which 
     follows chi square distribution with 4*nData-4 degrees of freedom.
  */
    if(ChiSumMin>ChistatSum) { 
      ChiSumMin=ChistatSum;
      df=4.0*(GV.nData)-4.0;
      lmfrqmin=lmfrq;
    }

    if(verboseflag==1) {
      fprintf(stdout,"%f  %f %f %d  %f \n",
	      ObsvHeader.fmax,lmfrq,ObsvHeader.FSmax,GV.nData,2.0*ChistatSum);
    }
  }

  LALFree(FaFbObsv);
  LALFree(FaFbTest);
  LALFree(GV.lmfrq);

  if(computeProbFlag) {
    prob = chi2cdf(df,2.0*ChiSumMin);
    if(prob<sigLevel) 
      rejection=1;
    fprintf(stdout,"%f  %f %f  %f %f %E %d\n",
	  ObsvHeader.fmax,lmfrqmin,ObsvHeader.FSmax,df,2.0*ChiSumMin,prob,rejection);
  } else {
    fprintf(stdout,"%f  %f %f  %f %f \n",
	    ObsvHeader.fmax,lmfrqmin,ObsvHeader.FSmax,df,2.0*ChiSumMin);
  }

  return 0;
}


/*-----------------------------------------------------------------*/
/*                                                                 */
/*-----------------------------------------------------------------*/

INT2 InitData(void)
{

  /* initialization */
  GV.nData=0;
  GV.indexOffset=0;
  GV.indexStepO=1;
  GV.indexStepT=1;
  GV.sFreq=0.0;
  GV.dFreq=0.0;
  GV.fmax=0.0;

  ObsvHeader.nData=0;
  ObsvHeader.sFreq=0.0;
  ObsvHeader.dFreq=0.0;
  ObsvHeader.fmax=0.0;

  TestHeader.nData=0;
  TestHeader.sFreq=0.0;  
  TestHeader.dFreq=0.0;
  TestHeader.fmax=0.0;

  return 0;
}


/*-----------------------------------------------------------------*/
/*                                                                 */
/*-----------------------------------------------------------------*/


INT2 HandleComArg(INT4 argc, CHAR *argv[]) 
{
  INT4 c, errflg = 0;
  
  /* scan through the list of arguments on the command line 
     and get the input data filename*/
  
  while (!errflg && ((c = getopt(argc, argv,"hpo:t:v:"))!=-1)) {
    switch (c) {
   case 'C':
      /* Verbose Output for debugging */
      computeProbFlag=atoi(optarg);
      break;
    case 'v':
      /* Verbose Output for debugging */
      verboseflag=atoi(optarg);
      break;
    case 'o':
      /* Name of observed data file */
      obsvdatafile=optarg;
      break;
    case 't':
      /* Name of test data file */
      testdatafile=optarg;
      break;
    case 'h':
      /* help */
      fprintf(stderr,"Usage: FstatShapeTest [-hv] [-o <>][-t <>]\n");
      fprintf(stderr,
	      "-o: <CHAR STRING:FaFbObsv.txt> File <filename> includes the observed data to be vetoed.\n");
      fprintf(stderr,
	      "-t: <CHAR STRING:FaFbTest.txt> File <filename> includes the veto signal: \n");
      fprintf(stderr,"FstatSphapeTest -v : Output in a verbose way. Mainly for debugging.\n");
      fprintf(stderr,"FstatSphapeTest -h : Show this help\n");
      fprintf(stderr,"Example: ./FstatShapeTest -o <ObservedDataFile> -t <VetoDataFile>\n");
      exit(0);
      break;
    default:
      /* unrecognized option */
      errflg++;
      fprintf(stderr,"Unrecognized option argument %c\n",c);
      exit(1);
      break;
    }
  }

return errflg;
}



/*-----------------------------------------------------------------*/
/*                                                                 */
/*-----------------------------------------------------------------*/

INT2 ReadHeader(void)
{ 
  INT4 errflg=0;
  INT4 nObsv,nTest;
  REAL8 sFreqO,dFreqO;  
  REAL8 sFreqT,dFreqT;
  REAL8 Aobsv,Bobsv,Cobsv;
  REAL8 Atest,Btest,Ctest;
  REAL8 fmaxObsv,fmaxTest;
  REAL8 Fmaxo,Fmaxt;
  const REAL8 errorTol=10e-7;
  FILE *fpobsv, *fptest;



  /* Counts the number of points in the cluster */
  /* read the hader of the observed data*/
  fpobsv = fopen(obsvdatafile,"r");
  if(fpobsv==NULL) {
    fprintf(stderr,"File open error in FstatShapeTest: %s\n",obsvdatafile);
    exit(1);
  }
  if(fscanf(fpobsv,"%d",&nObsv)==EOF) {    
    fprintf(stderr,"File format error in FstatShapeTest: %s\n",obsvdatafile);
    exit(1);
  }
  if(fscanf(fpobsv,"%lf %lf",&fmaxObsv,&Fmaxo)==EOF) {    
    fprintf(stderr,"File format error in FstatShapeTest: %s\n",obsvdatafile);
    exit(1);
  }
  if(fscanf(fpobsv,"%lf %lf",&sFreqO,&dFreqO)==EOF) {    
    fprintf(stderr,"File format error in FstatShapeTest: %s\n",obsvdatafile);
    exit(1);
  }
  if(fscanf(fpobsv,"%lf %lf %lf",&Aobsv,&Bobsv,&Cobsv)==EOF) {    
    fprintf(stderr,"File format error in FstatShapeTest: %s\n",obsvdatafile);
    exit(1);
  }
  fclose(fpobsv);



  /* read the hader of the test data*/
  fptest = fopen(testdatafile,"r");
  if(fptest==NULL) {
    fprintf(stderr,"File open error in FstatShapeTest: %s\n",testdatafile);
    exit(1);
  }
  if(fscanf(fptest,"%d",&nTest)==EOF) {    
    fprintf(stderr,"File format error in FstatShapeTest: %s\n",testdatafile);
    exit(1);
  }
  if(fscanf(fptest,"%lf %lf",&fmaxTest,&Fmaxt)==EOF) {    
    fprintf(stderr,"File format error in FstatShapeTest: %s\n",obsvdatafile);
    exit(1);
  }
  if(fscanf(fptest,"%lf %lf",&sFreqT,&dFreqT)==EOF) {    
    fprintf(stderr,"File format error in FstatShapeTest: %s\n",obsvdatafile);
    exit(1);
  }
  if(fscanf(fptest,"%lf %lf %lf",&Atest,&Btest,&Ctest)==EOF) {    
    fprintf(stderr,"File format error in FstatShapeTest: %s\n",obsvdatafile);
    exit(1);
  }
  fclose(fptest);



  if(nObsv>MAX_DATA_POINTS) {
    fprintf(stderr,
	    "Warning: Number of the data points exceeds MAX_DATA_POINTS.\n");
    exit(1);
  }
  if(nTest>MAX_DATA_POINTS) {
    fprintf(stderr,
	    "Warning: Number of the data points exceeds MAX_DATA_POINTS.\n");
    exit(1);
  }


  if(fabs(Aobsv-Atest)>fabs(Aobsv)*errorTol) { 
    fprintf(stderr,
	    "Warning: A of data is different from A of test.\n");
    return 1;
  }
  if(fabs(Bobsv-Btest)>fabs(Bobsv)*errorTol) {
    fprintf(stderr,
	    "Warning: A of data is different from B of test.\n");
    return 1;
  }
  if(fabs(Cobsv-Ctest)>fabs(Cobsv)*errorTol) { 
    fprintf(stderr,
	    "Warning: A of data is different from C of test.\n");
    return 1;
  }

  ObsvHeader.nData=nObsv;
  ObsvHeader.sFreq=sFreqO;
  ObsvHeader.dFreq=dFreqO;
  ObsvHeader.fmax=fmaxObsv;
  ObsvHeader.FSmax=Fmaxo;

  TestHeader.nData=nTest; 
  TestHeader.sFreq=sFreqT;
  TestHeader.dFreq=dFreqT;
  TestHeader.fmax=fmaxTest;
  TestHeader.FSmax=Fmaxt;

  GV.Acoeff=Aobsv;
  GV.Bcoeff=Bobsv;
  GV.Ccoeff=Cobsv;

  return errflg;
}

/*-----------------------------------------------------------------*/
                                                                    
/*-----------------------------------------------------------------*/
INT2 SummitFinder(void)
{
  REAL8 thr=0.5;
  INT4 iiter,jiter;
  INT4 nData;

  GV.lmfrq=(REAL8 *) LALCalloc(ObsvHeader.nData,sizeof(REAL8));

  nData=ObsvHeader.nData;
  jiter=0;
  for(iiter=0;iiter<nData;iiter++) {
    if(FaFbObsv[iiter].Fstat>=thr*ObsvHeader.FSmax) {
      GV.lmfrq[jiter]=FaFbObsv[iiter].freq;
      jiter++;
    }
  }
  GV.nSearch=jiter;

  return 0;
}


/*-----------------------------------------------------------------*/
                                                                    
/*-----------------------------------------------------------------*/
INT2 ShiftData(REAL8 lmfrq)
{
  /* This routine just shift the test data so that the frequency 
     where the maximum amplitude occurs coincide with that of the 
     data to be tested **/ 


  REAL8 sFreqT;
  REAL8 fmaxT;
  REAL8 freqoffset=0.0;

  /* pass the data to the local variables */
  fmaxT =TestHeader.fmax;
  sFreqT=TestHeader.sFreq;

  /* shift the data */
  freqoffset=lmfrq-fmaxT;
  TestHeader.sFreq=sFreqT+freqoffset;
  TestHeader.fmax=fmaxT+freqoffset;

  GV.fmax=ObsvHeader.fmax;


  return 0;
}


/*-----------------------------------------------------------------*/
  
/*-----------------------------------------------------------------*/
INT2 RearrangeData(void)
{
  REAL8 sFreqO,sFreqT,sFreq;
  REAL8 dFreqO,dFreqT,dFreqLCM;
  INT4   nObsv,nTest,nData;
  INT4 indexOffset=0;
  REAL8 errorTol=1.0*10e-5;

  /* pass the data to the local variables */
  nObsv =ObsvHeader.nData;
  sFreqO=ObsvHeader.sFreq;
  dFreqO=ObsvHeader.dFreq;

  nTest =TestHeader.nData;
  sFreqT=TestHeader.sFreq;
  dFreqT=TestHeader.dFreq;


  /* Maybe find LCM of dFreqT and dFreqO in the future */
  if(fabs(dFreqO-dFreqT)>dFreqO*errorTol) {
    fprintf(stderr,"The frequency resolution must be the same\n");
    return 1;
  }

  dFreqLCM=dFreqO;
  if(sFreqO-sFreqT!=errorTol) {
    if(fabs(fmod(sFreqO-sFreqT,dFreqLCM))>errorTol) {
      fprintf(stderr,"The starting frequency inconsistency\n");
      return 1;
    }
  }
  indexOffset= (INT4) myRound((sFreqO-sFreqT)/dFreqLCM);
  
  /* if dFreqLCM is really LCM, need to recompute nData */

  if(indexOffset>=0) {
    sFreq=sFreqO;
    if(nObsv+indexOffset>=nTest) {
      nData=nTest-indexOffset;
    } else {
      nData=nObsv;
    }
  } else {
    sFreq=sFreqT;
    if(nObsv+indexOffset>=nTest) {
      nData=nTest;
    } else {
      nData=nObsv+indexOffset;
    }
  }

  /* pass the data to the global variables */
  GV.nData=nData;
  GV.indexOffset=indexOffset;
  GV.sFreq=sFreq;
  GV.dFreq=dFreqLCM;

  return 0;
}

/*-----------------------------------------------------------------*/
/*                                                                 */
/*-----------------------------------------------------------------*/
INT4 myRound(REAL8 x)
{
  REAL8 sign=1.0;
  REAL8 roundedValue=0.0;
  REAL8 rmdr=0.0;

  if(x<0) sign=-1.0;
  roundedValue= floor(sign*x);
  rmdr=sign*x-roundedValue;
  if(rmdr>=0.5) 
    roundedValue=roundedValue+1.0;
  roundedValue=sign*roundedValue;

  return (INT4) roundedValue;
}



/*-----------------------------------------------------------------*/
/*                                                                 */
/*-----------------------------------------------------------------*/

INT2 ReadData(void)
{
  INT4 irec;
  CHAR buff[500];
  const INT4 Nheadlines=4;

  FILE *fpobsv,*fptest;

  fpobsv = fopen(obsvdatafile,"r");
  fptest = fopen(testdatafile,"r");

  /* Memory allocation */
  FaFbObsv=(FaFbdata *) LALCalloc(ObsvHeader.nData,sizeof(FaFbdata));
  FaFbTest=(FaFbdata *) LALCalloc(TestHeader.nData,sizeof(FaFbdata));


  /* skip the header */
  /* depend on data format specification */
  for(irec=0;irec<Nheadlines;irec++) {
    if ( fgets(buff,sizeof(buff),fpobsv) == NULL ) {
      fprintf (stderr, "\nfgets() failed!\n" );
      return 1;
    }
    if ( fgets(buff,sizeof(buff),fptest) == NULL ) {
      fprintf (stderr, "\nfgets() failed!\n" );
      return 1;
    }
  }

  /* data input begin */
  for(irec=0;irec<(ObsvHeader.nData);irec++) {
    int numread = fscanf(fpobsv,"%lf %lf %lf %lf %lf %lf",
	   &(FaFbObsv[irec].freq),
	   &(FaFbObsv[irec].RFa),
	   &(FaFbObsv[irec].IFa),
	   &(FaFbObsv[irec].RFb),
	   &(FaFbObsv[irec].IFb),
	   &(FaFbObsv[irec].Fstat)
	   );
    if ( numread != 6 ) {
      fprintf (stderr, "\nfscanf() failed to read 6 items from streadm 'fpobsv'\n");
      return 1;
    }
  }
  if(irec !=(ObsvHeader.nData)) {
    fprintf(stderr,"Data read error\n");
    return 1;
  }
  for(irec=0;irec<(TestHeader.nData);irec++) {
    int numread = fscanf(fptest,"%lf %lf %lf %lf %lf %lf",
	   &(FaFbTest[irec].freq),
	   &(FaFbTest[irec].RFa),
	   &(FaFbTest[irec].IFa),
	   &(FaFbTest[irec].RFb),
	   &(FaFbTest[irec].IFb),
	   &(FaFbTest[irec].Fstat)
	   );
    if ( numread != 6 ) {
      fprintf (stderr, "\nfscanf() failed to read 6 items from streadm 'fptest'\n");
      return 1;
    }

  }
  if(irec !=(TestHeader.nData)) {
    fprintf(stderr,"Data read error\n");
    return 1;
  }

  fclose(fpobsv);
  fclose(fptest);


  return 0;
}
/*-----------------------------------------------------------------*/
/*                                                                 */
/*-----------------------------------------------------------------*/

INT2 ConstructData(void)
{
  INT4 irec;
  INT4 jindO,jindT;
  INT4 indexOffset;
  INT4 indexStepO,indexStepT;


  /* Memory allocation */
  FaFbSubt=(FaFbdata *) LALCalloc(GV.nData,sizeof(FaFbdata));

  indexOffset=GV.indexOffset;
  indexStepO =GV.indexStepO;
  indexStepT =GV.indexStepT;

  if(indexOffset>=0) {
    for(irec=0;irec<(GV.nData);irec++) {
      jindO=irec*indexStepO;
      jindT=irec*indexStepT + indexOffset;

      FaFbSubt[irec].RFa = FaFbObsv[jindO].RFa-FaFbTest[jindT].RFa;
      FaFbSubt[irec].IFa = FaFbObsv[jindO].IFa-FaFbTest[jindT].IFa;
      FaFbSubt[irec].RFb = FaFbObsv[jindO].RFb-FaFbTest[jindT].RFb;
      FaFbSubt[irec].IFb = FaFbObsv[jindO].IFb-FaFbTest[jindT].IFb;

#ifdef DEBG_CONSTDATA 
      fprintf(stderr,"%f %f %f\n",
	      FaFbObsv[jindO].RFa,FaFbTest[jindT].RFa,
	      FaFbSubt[irec].RFa);
#endif
    }
  } else {
    for(irec=0;irec<(GV.nData);irec++) {
      jindO=irec*indexStepO - indexOffset;
      jindT=irec*indexStepT; 

      FaFbSubt[irec].RFa = FaFbObsv[jindO].RFa-FaFbTest[jindT].RFa;
      FaFbSubt[irec].IFa = FaFbObsv[jindO].IFa-FaFbTest[jindT].IFa;
      FaFbSubt[irec].RFb = FaFbObsv[jindO].RFb-FaFbTest[jindT].RFb;
      FaFbSubt[irec].IFb = FaFbObsv[jindO].IFb-FaFbTest[jindT].IFb;
#ifdef DEBG_CONSTDATA 
      fprintf(stderr,"%f %f %f\n",
	      FaFbObsv[jindO].RFa,FaFbTest[jindT].RFa,
	      FaFbSubt[irec].RFa);
#endif
    }
  }


  return 0;
}


/*-----------------------------------------------------------------*/
/*                                                                 */
/*-----------------------------------------------------------------*/


INT2 ComputeChi(void)
{ 
  INT4 irec;    /* counter */

  REAL8 Acoef,Bcoef,Ccoef,Dcoef;

  REAL8 RFa,IFa,RFb,IFb;
  REAL8 FaSq,FbSq,FaFb;
  REAL8 *Chistat;

  Acoef=GV.Acoeff;
  Bcoef=GV.Bcoeff;
  Ccoef=GV.Ccoeff;

  /* Memory allocation */
  Chistat  = (REAL8 *) LALCalloc(GV.nData,sizeof(REAL8));
  ChistatSum=0.0;

  for(irec=0;irec<(GV.nData);irec++) {


    /* reconstruct F */
    /* WARNING 1: 
       Here we assume that Fa and Fb are already 
       normaized by M (the number of sfts). 
       See the pulsargroup document or LALDemod document.

       This means that we use 
       F = (4.0/D)*(B*FaSq + A*FbSq - 2.0*C*FaFb); 
       instead of 
       F = (4.0/(M*D))*(B*FaSq + A*FbSq - 2.0*C*FaFb); 
    */

    /* WARNING 2: 
       Here we assume that the one sided noise power spectrum density 
       Sh is properly normalized. Namely, the Fa and Fb must be 
       multiplied by B*log(2.0) if one uses running median for an 
       estimate of Sh. B is the sample bias of the sample median 
       estimate.
    */
    
    RFa = FaFbSubt[irec].RFa;
    IFa = FaFbSubt[irec].IFa;
    RFb = FaFbSubt[irec].RFb;
    IFb = FaFbSubt[irec].IFb;
    

    Dcoef = Acoef*Bcoef - Ccoef*Ccoef;
    FaSq  = RFa*RFa + IFa*IFa;
    FbSq  = RFb*RFb + IFb*IFb;
    FaFb  = RFa*RFb + IFa*IFb;



    Chistat[irec] = (4.0/Dcoef)*(Bcoef*FaSq + Acoef*FbSq - 2.0*Ccoef*FaFb);
    ChistatSum += Chistat[irec];

    if(verboseflag==2)
      fprintf(stdout,"%f\n",Chistat[irec]);
  }

  LALFree(Chistat);
  LALFree(FaFbSubt);

  return 0;
}


/*-----------------------------------------------------------------*/
/*                                                                 */
/*-----------------------------------------------------------------*/

REAL8 chi2cdf(REAL8 df, REAL8 data)
{
  return gammaq(0.5*df,0.5*data);
}


/*-----------------------------------------------------------------*/
/*                                                                 */
/*-----------------------------------------------------------------*/


REAL8 gammaq(REAL8 a, REAL8 x)
{
  REAL8 gamser,gammcf,gln;

  if(x<0.0||a<=0.0) {
    fprintf(stderr,"Invalid arguments in gammaq\n");
    exit(1);
  }
  if(x<a+1.0) {
    gser(&gamser,a,x,&gln);
    return 1.0-gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return gammcf;
  }
}


void gser(REAL8 * gamser, REAL8 a, REAL8 x, REAL8 * gln)
{
  const INT4 ITMAX=100;
  const REAL8 EPS=3.0e-7;
  INT4 n;
  REAL8 sum,del,ap;

  *gln=gammln(a);
  if(x<=0.0) {
    if(x<0.0) {
      fprintf(stderr,"x less than 0 in gser");
      exit(1);
    }
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=1.0/a;
    sum=1.0/a;
    for(n=1;n<=ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if(fabs(del)<fabs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    fprintf(stderr,"a too large, ITMAX too small in gser");
    return;
  }
}


void gcf(REAL8 * gammcf, REAL8 a, REAL8 x, REAL8 * gln)
{
  const INT4 ITMAX=100;
  const REAL8 EPS=3.0e-7;
  const REAL8 FPMIN=1.0e-30;

  INT4 i;
  REAL8 an,b,c,d,del,h;

  *gln=gammln(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for(i=1;i<=ITMAX;i++) {
    an=-i*(i-a);
    b+=2.0;
    d=an*d+b;
    if(fabs(d)<FPMIN) d=FPMIN;
    c=b+an/c;
    if(fabs(c)<FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h*=del;
    if(fabs(del-1.0)<EPS) break;
  }
  if(i>ITMAX) {
    fprintf(stderr,"a too large, ITMAX too small in gcf");
    exit(1);
  }
  *gammcf=exp(-x+a*log(x)-(*gln))*h;
}



REAL8 gammln(REAL8 xx)
{
  REAL8 x,y,tmp,ser;
  static REAL8 cof[6]={76.18009172947146,
			-86.50532032941677,
			24.01409824083091,
			-1.231739572450155,
			0.1208650973866179e-2,
			-0.5395239384953e-5};
  INT4 j;

  y=xx;
  x=xx;
  tmp=x+5.5;
  tmp-=(x+0.5)*log(tmp);
  ser=1.000000000190015;
  for(j=0;j<=5;j++) ser+=cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}




