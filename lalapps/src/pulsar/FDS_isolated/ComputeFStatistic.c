/*********************************************************************************/
/*                    F-statistic generation code for known pulsars              */
/*                                                                               */
/*			      Y. Ioth, M.A. Papa, X. Siemens                     */
/*                                                                               */
/*                 Albert Einstein Institute/UWM - started September 2002        */
/*********************************************************************************/
#include <lal/LALDemod.h>

#include "ComputeFStatistic.h"
#include "rngmed.h"
#include "clusters.h"
#include "DopplerScan.h"

/* #define DEBG_FAFB                 */
/* #define DEBG_ESTSIGPAR */
/* #define DEBG_MAIN  */
/* #define DEBG_SGV */

/* If FILE_FILENAME is defined, then print the corresponding file  */
#define FILE_FSTATS
#define FILE_FMAX
/* #define FILE_FLINES */
/* #define FILE_FTXT */
/* #define FILE_PSD */
/* #define FILE_PSDLINES */
/* #define FILE_SPRNG */



FFT **SFTData=NULL;                 /* SFT Data for LALDemod */
DemodPar *DemodParams  = NULL;      /* Demodulation parameters for LALDemod */
LIGOTimeGPS *timestamps=NULL;       /* Time stamps from SFT data */
LALFstat Fstat;
INT4 lalDebugLevel=0,i,a,d,s,irec;
static LALStatus status;
AMCoeffs amc;
GlobalVariables GV;
REAL8 MeanOneOverSh=0.0;
REAL8 Alpha,Delta;
Clusters HFLines={0}, HPLines={0};
Clusters *highSpLines=&HPLines, *highFLines=&HFLines;
/* #ifdef FILE_FMAX     */
FILE *fpmax;
/* #endif */
/* #ifdef FILE_STATS */
FILE *fpstat;
/* #endif     */
REAL8 medianbias=1.0;

DopplerScanState_t thisScan;


int main(int argc,char *argv[]) 
{

  INT4 *maxIndex=NULL; /*  array that contains indexes of maximum of each cluster */
  DopplerPosition_t dopplerpos;
  DopplerScanInit_t scanInit;

  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;

  if (SetGlobalVariables(CommandLineArgs)) return 2;

  if (AllocateMem()) return 3;

  if (ReadSFTData()) return 4;

  /*  This fills-in highSpLines that are then used by NormaliseSFTRngMdn */
#if 0
  if (GV.noise!=1){
    if (EstimatePSDLines()) return 6;
  }
#endif

  if (NormaliseSFTDataRngMdn()) return 5;

#ifdef FILE_FMAX  
  /*   open file */
  if (!(fpmax=fopen("Fmax","w"))){
    fprintf(stderr,"in Main: unable to open Fmax file\n");
    return 2;
  }
#endif
#ifdef FILE_FSTATS  
  /*      open file */
  if (!(fpstat=fopen("Fstats","w"))){
    fprintf(stderr,"in Main: unable to open Fstats file\n");
    return 2;
  }
#endif
	
  scanInit.Alpha = GV.Alpha;
  scanInit.AlphaBand = GV.AlphaBand;
  scanInit.dAlpha = GV.dAlpha;
  scanInit.Delta = GV.Delta;
  scanInit.DeltaBand = GV.DeltaBand;
  scanInit.dDelta = GV.dDelta;

  /* use metric-grid? and if so, for what maximal mismatch? */
  scanInit.useMetric = GV.useMetric;
  scanInit.metricMismatch = GV.metricMismatch;
  scanInit.flipTiling = GV.flipTiling;

  scanInit.obsBegin = SFTData[0]->fft->epoch;
  scanInit.obsDuration = SFTData[GV.SFTno-1]->fft->epoch.gpsSeconds - scanInit.obsBegin.gpsSeconds + GV.tsft;
  /*   scanInit.obsDuration = 36000; */

  scanInit.fmax  = GV.Freq;
  if (GV.FreqBand > 0) scanInit.fmax += GV.FreqBand;
  scanInit.Detector = GV.Detector;



  InitDopplerScan ( &status, &thisScan, scanInit);
  if (status.statusCode != 0)
    {
      REPORTSTATUS( &status );
      return (-1);
    }
  
  
  while (1)
    {
      NextDopplerPos( &status, &dopplerpos, &thisScan );

      if (status.statusCode != 0)
	{
	  REPORTSTATUS( &status );
	  return (-1);
	}

      Alpha = dopplerpos.skypos.longitude;
      Delta = dopplerpos.skypos.latitude; 
      
      /* Have we scanned all DopplerPositions yet? */
      if (dopplerpos.finished)  
	break;
      
      if (CreateDemodParams()) return 6;
      /* loop over spin params */
      for(s=0;s<GV.SpinImax;s++)
	{
	  DemodParams->spinDwn[0]=GV.Spin+s*GV.dSpin;
	  LALDemod (&status, &Fstat, SFTData, DemodParams);
	  
	  /*  This fills-in highFLines that are then used by discardFLines */
	  if (GV.FreqImax > 5)
	    if (EstimateFLines()) return 6;
	  

#ifdef DEBG_MAIN
	  for(i=0;i < GV.FreqImax ;i++)
	    {
	      /* medianbias is 1 if GV.noise=1 */ 
	      fprintf(stdout,"%20.10f %e %20.17f %20.17f %20.17f\n",
		      GV.Freq+i*GV.dFreq,  DemodParams->spinDwn[0],
		      Alpha, Delta, 2.0*medianbias*Fstat.F[i]);
	    }
#endif
	  
	  /*  This fills-in highFLines  */
	  if (highFLines != NULL && highFLines->Nclusters > 0){
	    
	    maxIndex=(INT4 *)LALMalloc(highFLines->Nclusters*sizeof(INT4));
	    
	    /*  for every cluster writes the information about it in file Fstats */

	    if (writeFLines(maxIndex)){
	      fprintf(stderr, "%s: trouble making file Fstats\n", argv[0]);
	      return 6;
	    }
	  }
	  
	  if( GV.EstimSigParam &&(highFLines !=NULL) && (highFLines->Nclusters >0))
	    if(writeFaFb(maxIndex)) return 255;
	  
	  
	  if( GV.EstimSigParam &&(highFLines !=NULL) &&(highFLines->Nclusters >0)) {
	    if (EstimateSignalParameters(maxIndex)) return 7;
	  }
	  
	  if (PrintTopValues(/* thresh */ 0.0, /* max returned */ 1))
	    fprintf(stderr, "%s: trouble making files Fmax"
		    "and/or Fstats\n", argv[0]);	  
	  fflush(stderr);
	  
	  if (highFLines != NULL && highFLines->Nclusters > 0){
	    LALFree(maxIndex);
	  }

	  /* Set the number of the clusters detected to 0 at each iteration 
	     of the sky-direction and the spin down */
	  highFLines->Nclusters=0;

	} /* For GV.spinImax */
      
    } /*  while SkyPos */


#ifdef FILE_FMAX  
  fclose(fpmax);
#endif
#ifdef FILE_FSTATS  
  fclose(fpstat);
#endif
  if (Freemem()) return 7;

  return 0;
}


/*******************************************************************************/
/*  Note that there is a degeneracy where the shifts taken at the same time,  */
/*  psi -> psi+Pi/2 and Phi0 -> Phi0 + Pi/2,  */
/*  give the same A1,A2,A3,A4.  */
/*  */
/*******************************************************************************/


int EstimateSignalParameters(INT4 * maxIndex)
{
  INT4 irec,jrec;
  REAL8 A1,A2,A3,A4,A=amc.A,B=amc.B,C=amc.C,D=amc.D;
  REAL8 beta,Asq,detA,ampratio,A1test;
  REAL8 psi_mle,Phi0_mle,mu_mle;
  REAL8 h0mle,h0mleSq;
  REAL8 error_tol=1.0/pow(10,14);
  REAL8 norm;
  FILE * fpMLEParam;


#ifdef DEBG_ESTSIGPAR
  REAL8 Ftest;
  REAL8 A2test,A3test,A4test;
#endif


  if(!(fpMLEParam=fopen("ParamMLE.txt","w")))
     fprintf(stderr,"Error in EstimateSignalParameters: unable to open the file");


  norm=2.0*sqrt(GV.tsft)/(GV.tsft*GV.SFTno);


  for(jrec=0;jrec < highFLines->Nclusters ;jrec++)
    {

      irec=maxIndex[jrec];

      A1 =  2.0*( B * Fstat.Fa[irec].re - C * Fstat.Fb[irec].re) / D;
      A2 =  2.0*( A * Fstat.Fb[irec].re - C * Fstat.Fa[irec].re) / D;
      A3 = - 2.0*( B * Fstat.Fa[irec].im - C * Fstat.Fb[irec].im) / D;
      A4 = - 2.0*( A * Fstat.Fb[irec].im - C * Fstat.Fa[irec].im) / D;



      Asq = A1*A1 + A2*A2 + A3*A3 + A4*A4;
      detA = A1*A4-A2*A3;

      h0mle = 0.5*pow(
		      pow(((A1-A4)*(A1-A4)+(A2+A3)*(A2+A3)),0.25)+
		      pow(((A1+A4)*(A1+A4)+(A2-A3)*(A2-A3)),0.25)
		      ,2);

      h0mleSq = pow(h0mle,2.0);
      ampratio=Asq/h0mleSq;


      if(ampratio<0.25-error_tol||ampratio>2.0+error_tol) 
	{
	  fprintf(stderr,"Imaginary Cos[iota]; cannot compute parameters");
	  fprintf(stderr,"in the EstimateSignalParameters routine");
	  fprintf(stderr,"in ComputeFStatistic code");
	  fprintf(stderr,"Now exitting...");
	  /* 	  break; */
	  exit(1);
	}

      if(fabs(ampratio-0.25)<error_tol) {
	mu_mle =0.0;
      } else if(fabs(ampratio-2.0)<error_tol) {
	mu_mle = 1.0;
      } else {
	mu_mle = sqrt(-3.0+2.0*sqrt(2.0+ampratio));
      }

      if(detA<0) 
	mu_mle = - 1.0*mu_mle;


      if(Asq*Asq < 4.0*detA*detA)
	{
	  fprintf(stderr,"Imaginary beta; cannot compute parameters");
	  break;
	}

      /* Compute MLEs of psi and Phi0 up to sign of Cos[2*Phi0] */
      /* Make psi and Phi0 always in -Pi/2 to Pi/2 */ 
      beta  = (Asq + sqrt(Asq*Asq - 4.0*detA*detA))/(2.0*detA);
      psi_mle  = atan( (beta*A4-A1)/(beta*A3+A2) )/2.0;
      Phi0_mle  = atan( (A1-beta*A4)/(A3+beta*A2) )/2.0;


      A1test=h0mle*(0.5*(1+mu_mle*mu_mle)*cos(2.0*psi_mle)*cos(2.0*Phi0_mle)
		    -mu_mle*sin(2.0*psi_mle)*sin(2.0*Phi0_mle));

      /* Determine the sign of Cos[2*Phi0] */
      if(A1*A1test<0) {
	if(Phi0_mle>0) {
	  Phi0_mle=Phi0_mle - LAL_PI/2.0;
	} else {
	  Phi0_mle=Phi0_mle + LAL_PI/2.0;
	}
      }


#ifdef DEBG_ESTSIGPAR
      /* Reconstruct A1,A2,A3,A4. Compare them with the original values. */

      A1test=h0mle*(0.5*(1+mu_mle*mu_mle)*cos(2.0*psi_mle)*cos(2.0*Phi0_mle)
		    -mu_mle*sin(2.0*psi_mle)*sin(2.0*Phi0_mle));
      A2test=h0mle*(0.5*(1+mu_mle*mu_mle)*sin(2.0*psi_mle)*cos(2.0*Phi0_mle)
		    +mu_mle*cos(2.0*psi_mle)*sin(2.0*Phi0_mle));
      A3test=h0mle*(-0.5*(1+mu_mle*mu_mle)*cos(2.0*psi_mle)*sin(2.0*Phi0_mle)
		    -mu_mle*sin(2.0*psi_mle)*cos(2.0*Phi0_mle));
      A4test=h0mle*(-0.5*(1+mu_mle*mu_mle)*sin(2.0*psi_mle)*sin(2.0*Phi0_mle)
		    +mu_mle*cos(2.0*psi_mle)*cos(2.0*Phi0_mle));


      fprintf(stderr,"LALDemod_Estimate output: "
              "A1=%20.15f A2=%20.15f A3=%20.15f A4=%20.15f\n"
	      ,A1,A2,A3,A4);
      fprintf(stderr,"Reconstructed from MLE: "
              "A1=%20.15f A2=%20.15f A3=%20.15f A4=%20.15f !!!!\n\n",
	      A1test,A2test,A3test,A4test);
      fflush(stderr);


      if(fabs(A1-A1test)>fabs(A1)/(10e5)){ 
	fprintf(stderr,"Something is wrong with Estimate A1\n");
	fprintf(stderr,"Frequency index %d, %lf (Hz),A1=%f,A1test=%f\n",
		irec,GV.Freq+irec*GV.dFreq,A1,A1test);
	fprintf(stderr,"relative error Abs((A1-A1test)/A1)=%lf\n",
		fabs(A1-A1test)/fabs(A1));
	exit(1);
      }
      if(fabs(A2-A2test)>fabs(A2)/(10e5)){ 
	fprintf(stderr,"Something is wrong with Estimate A2\n");
	fprintf(stderr,"Frequency index %d, %lf (Hz),A2=%f,A2test=%f\n",
		irec,GV.Freq+irec*GV.dFreq,A2,A2test);
	fprintf(stderr,"relative error Abs((A2-A2test)/A2)=%lf\n",
		fabs(A2-A2test)/fabs(A2));
	exit(1);
      }
      if(fabs(A3-A3test)>fabs(A3)/(10e5)){ 
	fprintf(stderr,"Something is wrong with Estimate A3\n");
	fprintf(stderr,"Frequency index %d, %lf (Hz),A3=%f,A3test=%f\n",
		irec,GV.Freq+irec*GV.dFreq,A3,A3test);
	fprintf(stderr,"relative error Abs((A3-A3test)/A3)=%lf\n",
		fabs(A3-A3test)/fabs(A3));
	exit(1);
      }
      if(fabs(A4-A4test)>fabs(A4)/(10e5)){ 
	fprintf(stderr,"Something is wrong with Estimate A4\n");
	fprintf(stderr,"Frequency index %d, %lf (Hz),A4=%f,A4test=%f\n",
		irec,GV.Freq+irec*GV.dFreq,A1,A1test);
	fprintf(stderr,"relative error Abs((A4-A4test)/A4)=%lf\n",
		fabs(A4-A4test)/fabs(A4));
	exit(1);
      }

      
      /* Reconstruct F. Compare it with the original value. */
      Ftest=(A*(A1*A1+A3*A3)+B*(A2*A2+A4*A4)+2.0*C*(A1*A2+A3*A4))/
	4.0*4.0/GV.SFTno;


      if(fabs(Fstat.F[irec] - Ftest)> fabs(Ftest)/10e5){ 
	fprintf(stderr,"Something is wrong with Estimate in F\n");
	fprintf(stderr,"Frequency index %d, %lf (Hz),F=%f,Ftest=%f\n",
		irec,GV.Freq+irec*GV.dFreq,Fstat.F[irec],Ftest);
	fprintf(stderr,"relative error Abs((F-Ftest)/Ftest)=%lf\n",
		fabs(Fstat.F[irec]-Ftest)/fabs(Ftest));
	exit(1);
      }
#endif


      /* normalization */
      h0mle=h0mle*norm;


      /* For the real data, we need to multiply long(2.0) */
      /* Because we use running median to estimate the S_h. */
      /* if(GV.noise!=1) 
	h0mle=h0mle*sqrt(medianbias);
      */
      /* medianbias is 1 when GV.noise==1 */
      h0mle=h0mle*sqrt(medianbias);

#ifdef DEBG_ESTSIGPAR
      {double hp,hc,ds;
      hp=(1.0+mu_mle*mu_mle)*h0mle/2.0;
      hc=mu_mle*h0mle;
      ds=GV.SFTno*GV.tsft/2.0*(hp*hp*((A+B)/2.0+(A-B)/2.0*cos(4.0*psi_mle)
				   +C*sin(4.0*psi_mle))
			    +hc*hc*((A+B)/2.0-(A-B)/2.0*cos(4.0*psi_mle)
				    -C*sin(4.0*psi_mle)));
      fprintf(stderr,"A=%f,B=%f,C=%f,f=%f,h0=%f,F=%f\n",
	      A,B,C,GV.Freq+irec*GV.dFreq,h0mle,Fstat.F[irec]*medianbias);
      }
#endif

      /* Note that we print out MLE of 2.0*Phi0_JKS */
      /* because Phi0_PULGROUPDOC=2.0*Phi0_JKS */
      /* and Phi0_PULGROUPDOC is the one used in In.data. */
 
      /* medianbias is 1 if GV.noise==1 */
      fprintf(fpMLEParam,"%16.8lf %22E", GV.Freq+irec*GV.dFreq,2.0*medianbias*Fstat.F[irec]);


      fprintf(fpMLEParam,"  %10.6f",(1.0+mu_mle*mu_mle)*h0mle/2.0);
      fprintf(fpMLEParam,"  %10.6f",mu_mle*h0mle);
      fprintf(fpMLEParam,"  %10.6f",psi_mle);
      fprintf(fpMLEParam,"  %10.6f",2.0*Phi0_mle);
      fprintf(fpMLEParam,"\n");
    }

  fclose(fpMLEParam);

  return 0;

} /* EstimateSignalParameters() */


/*******************************************************************************/

/* Write the Fa and Fb for the later use of Fstatistic Shape test */
/* the explicit format specifier like %22.12f looks ugly. */
int writeFaFb(INT4 *maxIndex)
{
  INT4 irec,jrec;
  INT4 index,krec=0;
  CHAR filebasename[]="FaFb"; /* Base of the output file name */
  CHAR filename[256];         /* Base of the output file name */
  CHAR noiseswitch[16];
  CHAR clusterno[16];
  INT4 N;
  FILE * fp;
  REAL8 bias=1.0;

  sprintf(noiseswitch,"%02d",GV.noise);
  strcat(filebasename,noiseswitch);

  /*
  if(GV.noise!=1) 
    bias=sqrt(medianbias);
  */
  /* medianbias=1 if GV.noise==1 */
  bias=sqrt(medianbias);


  for (irec=0;irec<highFLines->Nclusters;irec++){
    sprintf(clusterno,".%03d",irec+1);
    strcpy(filename,filebasename);
    strcat(filename,clusterno);

    if((fp=fopen(filename,"w"))==NULL) {
      fprintf(stderr,"Unable to open a file %s\n",filename);
      return 1;
    }

    N=highFLines->NclustPoints[irec];


    /* The header contains */
    /* N the number of points in the cluster */
    /* the frequency where the maximum amplitude is */
    /* A,B,C coefficients */
    index=highFLines->Iclust[krec];

    fprintf(fp,"%10d\n",N);
    fprintf(fp,"%22.12f %22.12f\n",
	    GV.Freq+maxIndex[irec]*GV.dFreq,
	    Fstat.F[maxIndex[irec]]*bias*bias);
    fprintf(fp,"%22.12f %22.12f\n",GV.Freq+index*GV.dFreq,GV.dFreq);
    fprintf(fp,"%22.12f %22.12f %22.12f\n",amc.A,amc.B,amc.C);



    /* WARNING 1: 
       Here we make it sure that Fa and Fb are already 
       normaized by M (the number of sfts). 
       See the pulsargroup document or LALDemod document.

       This means that we use in the FstatShapeTest code  
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


    /* WARNING 3:
       The information stored in the FaFb file will be used by 
       FstatShapeTest code. The format of the file must not be 
       changed without the appropriate change in FstatShapeTest 
       code. 
     */ 

    for(jrec=0;jrec<N;jrec++) {
      index=highFLines->Iclust[krec];
      krec++;

#ifdef DEBG_FAFB
      fprintf(fp,"%22.16f %22.16f "
                 "%E %20.17f %20.17f "
                 "%22.16f %22.16f %22.16f %22.16f %22.16f %22.16f %22.16f\n",
	      GV.Freq+index*GV.dFreq,Fstat.F[index]*bias*bias,
	      DemodParams->spinDwn[0], Alpha, Delta,
	      Fstat.Fa[index].re/sqrt(GV.SFTno)*bias,
	      Fstat.Fa[index].im/sqrt(GV.SFTno)*bias,
	      Fstat.Fb[index].re/sqrt(GV.SFTno)*bias,
	      Fstat.Fb[index].im/sqrt(GV.SFTno)*bias,
	      amc.A,amc.B,amc.C);
#else
      /* Freqency, Re[Fa],Im[Fa],Re[Fb],Im[Fb], F */
      fprintf(fp,"%22.16f %22.12f %22.12f %22.12f %22.12f %22.12f\n",
	      GV.Freq+index*GV.dFreq,
	      Fstat.Fa[index].re/sqrt(GV.SFTno)*bias,
	      Fstat.Fa[index].im/sqrt(GV.SFTno)*bias,
	      Fstat.Fb[index].re/sqrt(GV.SFTno)*bias,
	      Fstat.Fb[index].im/sqrt(GV.SFTno)*bias,
	      Fstat.F[index]*bias*bias);
#endif


    }
    fclose(fp);
  }

  return 0;
}



/*******************************************************************************/

int CreateDemodParams(void)
{
  CSParams *csParams  = NULL;        /* ComputeSky parameters */
  BarycenterInput baryinput;         /* Stores detector location and other barycentering data */
  EphemerisData *edat=NULL;          /* Stores earth/sun ephemeris data for barycentering */
  EarthState earth;
  EmissionTime emit;
  AMCoeffsParams *amParams;
  LIGOTimeGPS *midTS=NULL;           /* Time stamps for amplitude modulation coefficients */
  LALLeapSecFormatAndAcc formatAndAcc = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
  INT4 leap;

  INT4 k;
  
  edat=(EphemerisData *)LALMalloc(sizeof(EphemerisData));
  (*edat).ephiles.earthEphemeris = GV.EphemEarth;     
  (*edat).ephiles.sunEphemeris = GV.EphemSun;         

  LALLeapSecs(&status,&leap,&timestamps[0],&formatAndAcc);
  (*edat).leap=leap;

  LALInitBarycenter(&status, edat);               /* Reads in ephemeris files */
 
/* Detector location: MAKE INTO INPUT!!!!! */
  baryinput.site.location[0]=GV.Detector.location[0]/LAL_C_SI;
  baryinput.site.location[1]=GV.Detector.location[1]/LAL_C_SI;
  baryinput.site.location[2]=GV.Detector.location[2]/LAL_C_SI;
  baryinput.alpha=Alpha;
  baryinput.delta=Delta;
  baryinput.dInv=0.e0;

/* amParams structure to compute a(t) and b(t) */

/* Allocate space for amParams stucture */
/* Here, amParams->das is the Detector and Source info */
  amParams = (AMCoeffsParams *)LALMalloc(sizeof(AMCoeffsParams));
  amParams->das = (LALDetAndSource *)LALMalloc(sizeof(LALDetAndSource));
  amParams->das->pSource = (LALSource *)LALMalloc(sizeof(LALSource));
/* Fill up AMCoeffsParams structure */
  amParams->baryinput = &baryinput;
  amParams->earth = &earth; 
  amParams->edat = edat;
  amParams->das->pDetector = &GV.Detector; 
  amParams->das->pSource->equatorialCoords.latitude = Delta;
  amParams->das->pSource->equatorialCoords.longitude = Alpha;
  amParams->das->pSource->orientation = 0.0;
  amParams->das->pSource->equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  amParams->polAngle = amParams->das->pSource->orientation ; /* These two have to be the same!!!!!!!!!*/
  amParams->leapAcc=formatAndAcc.accuracy;

 /* Mid point of each SFT */
   midTS = (LIGOTimeGPS *)LALCalloc(GV.SFTno,sizeof(LIGOTimeGPS));
   for(k=0; k<GV.SFTno; k++)
     { 
       REAL8 teemp=0.0;
       LALGPStoFloat(&status,&teemp, &(timestamps[k]));
       teemp += 0.5*GV.tsft;
       LALFloatToGPS(&status,&(midTS[k]), &teemp);
     }
   
   LALComputeAM(&status, &amc, midTS, amParams); 

/* ComputeSky stuff*/
  csParams=(CSParams *)LALMalloc(sizeof(CSParams));
  csParams->tGPS=timestamps;  
  csParams->skyPos=(REAL8 *)LALMalloc(2*sizeof(REAL8));
  csParams->mObsSFT=GV.SFTno;     /* Changed this from GV.mobssft !!!!!! */
  csParams->tSFT=GV.tsft;
  csParams->edat=edat;
  csParams->baryinput=&baryinput;
  csParams->spinDwnOrder=1;
  csParams->skyPos[0]=Alpha;
  csParams->skyPos[1]=Delta;
  csParams->earth = &earth;
  csParams->emit = &emit;

/* Finally, DemodParams */
  DemodParams->amcoe=&amc;
  DemodParams->spinDwnOrder=1;
  DemodParams->SFTno=GV.SFTno;

  DemodParams->f0=GV.Freq;
  DemodParams->imax=GV.FreqImax;
  DemodParams->df=GV.dFreq;

  DemodParams->Dterms=GV.Dterms;
  DemodParams->ifmin=GV.ifmin;

  DemodParams->returnFaFb = GV.EstimSigParam;

  ComputeSky(&status,DemodParams->skyConst,0,csParams);        /* compute the */
						               /* "sky-constants" A and B */
  LALFree(midTS);

  LALFree(csParams->skyPos);
  LALFree(csParams);

  LALFree(amParams->das->pSource);
  LALFree(amParams->das);
  LALFree(amParams);

  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);

  return 0;
}

/*******************************************************************************/
int AllocateMem(void)
{
  INT4 k;


  /* Allocate space for AMCoeffs */
  amc.a = NULL;
  amc.b = NULL;
  LALSCreateVector(&status, &(amc.a), (UINT4) GV.SFTno);
  LALSCreateVector(&status, &(amc.b), (UINT4) GV.SFTno);

  /* Allocate DemodParams structure */
  DemodParams=(DemodPar *)LALMalloc(sizeof(DemodPar));
  
  /* space for sky constants */
  /* Based on maximum index for array of as and bs sky constants as from ComputeSky.c */
  k=4*(GV.SFTno-1)+4; 
  DemodParams->skyConst=(REAL8 *)LALMalloc(k*sizeof(REAL8));

  /* space for spin down params */
  DemodParams->spinDwn=(REAL8 *)LALMalloc(sizeof(REAL8));
  
  return 0;
}

/*******************************************************************************/

/*  for every cluster writes the information about it in file Fstats */
/*  precisely it writes: */
/*  fr_max alpha delta N_points_of_cluster mean std max (of 2F) */
int writeFLines(INT4 *maxIndex){

  INT4 i,j,j1,j2,k,N;
  REAL8 max,log2,mean,var,std,R,fr;
  INT4 imax;
  INT4 err;

  log2=medianbias;
 
  j1=0;
  j2=0;

  for (i=0;i<highFLines->Nclusters;i++){
    N=highFLines->NclustPoints[i];
    
    /*  determine maximum of j-th cluster */
    /*  and compute mean */
    max=0.0;
    imax=0;
    mean=0.0;
    std=0.0;
    for (j=0;j<N;j++){
      R=2.0*log2*highFLines->clusters[j1];
      k=highFLines->Iclust[j1];
      j1=j1+1;
      mean=mean+R;
      if( R > max){
	max=R;
	imax=k;
      }
    }/*  end j loop over points of i-th cluster  */
    /*  and start again to compute variance */
    maxIndex[i]=imax;
    mean=mean/N;
    var=0.0;
    for (j=0;j<N;j++){
      R=2.0*log2*highFLines->clusters[j2];
      j2=j2+1;
      var=var+(R-mean)*(R-mean);
    }/*  end j loop over points of i-th cluster  */
    var=var/N;
    std=sqrt(var);
    fr=GV.Freq+imax*GV.dFreq;
#ifdef FILE_FSTATS  
/*    print the output */
  err=fprintf(fpstat,"%16.12f %10.8f %10.8f    %d %10.5f %10.5f %10.5f\n",fr,
	      Alpha, Delta, N, mean, std, max);

  if (err<=0) {
    fprintf(stderr,"writeFLines couldn't print to Fstas!\n");
    return 4;
  }
#endif

  }/*  end i loop over different clusters */

  return 0;
}

/*******************************************************************************/

int NormaliseSFTData(void)
{
  INT4 k,j;                         /* loop indices */
  INT4 nbins=GV.ifmax-GV.ifmin+1;   /* Number of points in SFT's */
  REAL8 SFTsqav;                  /* Average of Square of SFT */
  REAL8 B;                        /* SFT Bandwidth */
  REAL8 deltaT,N;


  /* loop over each SFTs */
  for (k=0;k<GV.SFTno;k++)         
    {
      
      SFTsqav=0.0;
      /* loop over SFT data to estimate noise */
      for (j=0;j<nbins;j++)               
	  {
	    SFTsqav=SFTsqav+
	      SFTData[k]->fft->data->data[j].re * SFTData[k]->fft->data->data[j].re+
	      SFTData[k]->fft->data->data[j].im * SFTData[k]->fft->data->data[j].im;
	  }
      SFTsqav=SFTsqav/(1.0*nbins);              /* Actual average of Square of SFT */
      MeanOneOverSh=2.0*GV.nsamples*GV.nsamples/(SFTsqav*GV.tsft)+MeanOneOverSh;      

      N=1.0/sqrt(2.0*(REAL8)SFTsqav);

      /* signal only case */  
      if(GV.noise == 1)
	{
	  B=(1.0*GV.nsamples)/(1.0*GV.tsft);
	  deltaT=1.0/(2.0*B);
	  N=deltaT/sqrt(GV.tsft);
	}

      /* loop over SFT data to Normalise it*/
      for (j=0;j<nbins;j++)               
	{
	  SFTData[k]->fft->data->data[j].re = N*SFTData[k]->fft->data->data[j].re; 
	  SFTData[k]->fft->data->data[j].im = N*SFTData[k]->fft->data->data[j].im;
	}
    } 
	
   MeanOneOverSh=MeanOneOverSh/(1.0*GV.SFTno); 
  return 0;
}

/*******************************************************************************/

int ReadSFTData(void)
{
  INT4 fileno=0,ndeltaf,offset;
  FILE *fp;
  size_t errorcode;

  SFTData=(FFT **)LALMalloc(GV.SFTno*sizeof(FFT *));
  timestamps=(LIGOTimeGPS *)LALMalloc(GV.SFTno*sizeof(LIGOTimeGPS));

  for (fileno=0;fileno<GV.SFTno;fileno++)
    {
      /* open FIRST file and get info from it*/

      fp=fopen(GV.filelist[fileno],"r");
      if (fp==NULL) 
	{
	  fprintf(stderr,"Weird... %s doesn't exist!\n",GV.filelist[fileno]);
	  return 1;
	}
      /* Read in the header from the file */
      errorcode=fread((void*)&header,sizeof(header),1,fp);
      if (errorcode!=1) 
	{
	  fprintf(stderr,"No header in data file %s\n",GV.filelist[fileno]);
	  return 1;
	}

      /* Check that data is correct endian order */
      if (header.endian!=1.0)
	{
	  fprintf(stderr,"First object in file %s is not (double)1.0!\n",GV.filelist[fileno]);
	  fprintf(stderr,"It could be a file format error (big/little\n");
	  fprintf(stderr,"endian) or the file might be corrupted\n\n");
	  return 2;
	}
    
      /* Check that the time base is positive */
      if (header.tbase<=0.0)
	{
	  fprintf(stderr,"Timebase %f from data file %s non-positive!\n",
		  header.tbase,GV.filelist[fileno]);
	  return 3;
	}
        
      /* Check that are frequency bins needed are in data set */
      if (GV.ifmin<header.firstfreqindex || 
	  GV.ifmax>header.firstfreqindex+header.nsamples) 
	{
	  fprintf(stderr,"Freq index range %d->%d not in %d to %d (file %s)\n",
		GV.ifmin,GV.ifmax,header.firstfreqindex,
		  header.firstfreqindex+header.nsamples,GV.filelist[fileno]);
	  return 4;
	}
      /* Put time stamps from file into array */
      timestamps[fileno].gpsSeconds=header.gps_sec;
      timestamps[fileno].gpsNanoSeconds=header.gps_nsec;

      /* Move forward in file */
      offset=(GV.ifmin-header.firstfreqindex)*2*sizeof(REAL4);
      errorcode=fseek(fp,offset,SEEK_CUR);
      if (errorcode) 
	{
	  perror(GV.filelist[fileno]);
	  fprintf(stderr,"Can't get to offset %d in file %s\n",offset,GV.filelist[fileno]);
	  return 5;
	}

      /* Make data structures */
      ndeltaf=GV.ifmax-GV.ifmin+1;
      SFTData[fileno]=(FFT *)LALMalloc(sizeof(FFT));
      SFTData[fileno]->fft=(COMPLEX8FrequencySeries *)LALMalloc(sizeof(COMPLEX8FrequencySeries));
      SFTData[fileno]->fft->data=(COMPLEX8Vector *)LALMalloc(sizeof(COMPLEX8Vector));
      SFTData[fileno]->fft->data->data=(COMPLEX8 *)LALMalloc(ndeltaf*sizeof(COMPLEX8));

      /* Fill in actual SFT data, and housekeeping */
      errorcode=fread((void*)(SFTData[fileno]->fft->data->data), sizeof(COMPLEX8), ndeltaf, fp);
      if (errorcode!=ndeltaf){
	perror(GV.filelist[fileno]);
	fprintf(stderr, "The SFT data was truncated.  Only read %d not %d complex floats\n", errorcode, ndeltaf);
	return 6;
      }
      SFTData[fileno]->fft->epoch=timestamps[fileno];
      SFTData[fileno]->fft->f0=GV.ifmin*GV.df;
      SFTData[fileno]->fft->deltaF=GV.df;
      SFTData[fileno]->fft->data->length=ndeltaf;

      fclose(fp);     /* close file */
    
    }
  return 0;  
}

/*******************************************************************************/

INT4 SetGlobalVariables(struct CommandLineArgsTag CLA)
{

  CHAR command[256];
  FILE *fp;
  size_t errorcode;
  REAL8 df;                         /* freq resolution */
  INT4 fileno=0;   
  glob_t globbuf;


  strcpy(GV.EphemEarth,CLA.EphemDir);
  strcat(GV.EphemEarth,"/earth");
  strcat(GV.EphemEarth,CLA.EphemYear);
  strcat(GV.EphemEarth,".dat");

  strcpy(GV.EphemSun,CLA.EphemDir);
  strcat(GV.EphemSun,"/sun");
  strcat(GV.EphemSun,CLA.EphemYear);
  strcat(GV.EphemSun,".dat");

  /* *** Make sure the e-files are really there *** */
      fp=fopen(GV.EphemEarth,"r");
      if (fp==NULL) 
	{
	  fprintf(stderr,"Could not find %s\n",GV.EphemEarth);
	  return 1;
	}
      fclose(fp);
      fp=fopen(GV.EphemSun,"r");
      if (fp==NULL) 
	{
	  fprintf(stderr,"Could not find %s\n",GV.EphemSun);
	  return 1;
	}
      fclose(fp);
  /* ********************************************** */

  strcpy(command,CLA.DataDir);
  strcat(command,"/*");
  strcat(command,CLA.BaseName);
  strcat(command,"*");
  
  globbuf.gl_offs = 1;
  glob(command, GLOB_ERR, NULL, &globbuf);

  /* read file names -- MUST NOT FORGET TO PUT ERROR CHECKING IN HERE !!!! */

  if(globbuf.gl_pathc==0)
    {
      fprintf(stderr,"No SFTs in directory %s ... Exiting.\n",CLA.DataDir);
      return 1;
    }

  while ((UINT4)fileno < globbuf.gl_pathc) 
    {
      strcpy(GV.filelist[fileno],globbuf.gl_pathv[fileno]);
      fileno++;
      if (fileno > MAXFILES)
	{
	  fprintf(stderr,"Too many files in directory! Exiting... \n");
	  return 1;
	}
    }
  globfree(&globbuf);

  GV.SFTno=fileno; /* remember this is 1 more than the index value */

  GV.IFO=CLA.IFO;

  /* initialize detector */
  if(GV.IFO == 0) GV.Detector=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if(GV.IFO == 1) GV.Detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
  if(GV.IFO == 2) GV.Detector=lalCachedDetectors[LALDetectorIndexLHODIFF];
  if(GV.IFO == 3)
    {
        if (CreateDetector(&(GV.Detector))) return 5;
    }


  /* open FIRST file and get info from it*/
  fp=fopen(GV.filelist[0],"r");
  /* read in the header from the file */
  errorcode=fread((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1) 
    {
      fprintf(stderr,"No header in data file %s\n",GV.filelist[0]);
      return 1;
    }

  /* check that data is correct endian order */
  if (header.endian!=1.0)
    {
      fprintf(stderr,"First object in file %s is not (double)1.0!\n",GV.filelist[0]);
      fprintf(stderr,"It could be a file format error (big/little\n");
      fprintf(stderr,"endian) or the file might be corrupted\n\n");
      return 2;
    }
  fclose(fp);

  GV.Ti=header.gps_sec;  /* INITIAL TIME */

  /* open LAST file and get info from it*/
  fp=fopen(GV.filelist[fileno-1],"r");
  /* read in the header from the file */
  errorcode=fread((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1) 
    {
      fprintf(stderr,"No header in data file %s\n",GV.filelist[fileno-1]);
      return 1;
    }
  /* check that data is correct endian order */
  if (header.endian!=1.0)
    {
      fprintf(stderr,"First object in file %s is not (double)1.0!\n",GV.filelist[fileno-1]);
      fprintf(stderr,"It could be a file format error (big/little\n");
      fprintf(stderr,"endian) or the file might be corrupted\n\n");
      return 2;
    }
  fclose(fp);

  GV.Tf=header.gps_sec+header.tbase;  /* FINAL TIME */

  GV.tsft=header.tbase;  /* Time baseline of SFTs */
    
  /* variables for starting demodulation frequency, band and resolution */
  GV.Freq=CLA.Freq;
  GV.FreqBand=CLA.FreqBand;
  GV.dFreq=CLA.dFreq;
  /* if user has not input demodulation frequency resolution; set to 1/Tobs */
  if( GV.dFreq == 0.0 ) GV.dFreq=1.0/(2.0*header.tbase*GV.SFTno);
  GV.FreqImax=(int)(GV.FreqBand/GV.dFreq+.5)+1;  /*Number of frequency values to calculate F for */
    
  /* variables for starting spindown, band and resolution */
  GV.Spin=CLA.Spin;
  GV.SpinBand=CLA.SpinBand;
  GV.dSpin=CLA.dSpin;
  /* if user has not input demodulation frequency resolution; set to 1/Tobs */
  if( GV.dSpin == 0.0 ) GV.dSpin=1.0/(2.0*header.tbase*GV.SFTno*(GV.Tf-GV.Ti));

  GV.SpinImax=(int)(GV.SpinBand/GV.dSpin+.5)+1;  /*Number of frequency values to calculate F for */

  /* sky position, band and resolution */
  GV.Alpha=CLA.Alpha;
  GV.AlphaBand=CLA.AlphaBand;
  GV.dAlpha=CLA.dAlpha;

  GV.Delta=CLA.Delta;
  GV.DeltaBand=CLA.DeltaBand;
  GV.dDelta=CLA.dDelta;

  GV.nsamples=header.nsamples;    /* # of freq. bins */

  /* frequency resolution: used only for printing! */
  df=(1.0)/(1.0*header.tbase);
  GV.df=df;

  GV.ifmax=ceil((1.0+DOPPLERMAX)*(GV.Freq+GV.FreqBand)*GV.tsft)+CLA.Dterms;
  GV.ifmin=floor((1.0-DOPPLERMAX)*GV.Freq*GV.tsft)-CLA.Dterms;

  GV.Dterms=CLA.Dterms;

  GV.noise=CLA.noise;

  GV.EstimSigParam=CLA.EstimSigParam;
  GV.Fthreshold=CLA.Fthreshold;

  /* allocate F-statistic arrays */
  Fstat.F =(REAL8*)LALMalloc(GV.FreqImax*sizeof(REAL8));
  if(GV.EstimSigParam) 
    {
      Fstat.Fa =(COMPLEX16*)LALMalloc(GV.FreqImax*sizeof(COMPLEX16));
      Fstat.Fb =(COMPLEX16*)LALMalloc(GV.FreqImax*sizeof(COMPLEX16));
    } else {
      Fstat.Fa = NULL;
      Fstat.Fb = NULL;
    }

  /* Tell the user what we have arrived at */
#ifdef DEBG_SGV
    fprintf(stdout,"\n");
    fprintf(stdout,"# SFT time baseline:                  %f min\n",header.tbase/60.0);
    fprintf(stdout,"# SFT freq resolution:                %f Hz\n",df);
    fprintf(stdout,"# Starting search frequency:          %f Hz\n",GV.Freq);
    fprintf(stdout,"# Demodulation frequency band:        %f Hz\n",GV.FreqBand);
    fprintf(stdout,"# no of SFT in a DeFT:                %f\n",ceil((1.0*(GV.Tf - GV.Ti))/header.tbase));
    fprintf(stdout,"# Actual # of SFTs:                   %d\n",GV.SFTno);
    fprintf(stdout,"# ==> DeFT baseline:                  %f hours\n",(GV.Tf - GV.Ti)/3600.0);
#endif

  return 0;  
}

/*******************************************************************************/

INT4 ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA) 
{
  INT4 c, errflg = 0;
  INT4 option_index = 0;

  /* ---------------------------------------------------------------------- */
  /* Initialize default values */
  CLA->Dterms=16;

  CLA->Freq=0.0;
  CLA->dFreq=0.0;
  CLA->FreqBand=0.0;

  CLA->Alpha=0.0;
  CLA->AlphaBand=0.0;
  CLA->Delta=0.0;
  CLA->DeltaBand=0.0;
  CLA->Spin=0.0;            
  CLA->SpinBand=0.0;

  /* these will not be used if the metric approach will be used: cmd-line option --use-mismatch -M */
  CLA->dDelta=0.001;
  CLA->dAlpha=0.001;
  CLA->dSpin=0.0;

  CLA->DataDir="";
  CLA->EphemYear="";
  CLA->EphemDir="";
  CLA->IFO=-1;
  CLA->noise=0;

  CLA->EstimSigParam=0;
  CLA->Fthreshold=10.0;
  CLA->BaseName="SFT";

  GV.useMetric = DOPPLER_MANUAL;
  GV.metricMismatch = 0.02;
  GV.flipTiling = 0;
  /* ---------------------------------------------------------------------- */

  while (1)
    {
      static struct option long_options[] =
	{
	  {"Freq", 		required_argument, 0, 	'f'},
	  {"FreqBand", 		required_argument, 0, 	'b'},
	  {"dFreq", 		required_argument,0, 	'r'},
	  {"Delta", 		required_argument, 0, 	'd'},
	  {"dDelta", 		required_argument, 0, 	'g'},
	  {"DeltaBand", 	required_argument, 0, 	'c'},
	  {"Alpha", 		required_argument, 0, 	'a'},
	  {"dAlpha", 		required_argument, 0, 	'l'},
	  {"AlphaBand", 	required_argument, 0, 	'z'},
	  {"Dterms", 		required_argument, 0, 	't'},
	  {"EphemYear", 	required_argument, 0, 	'y'},
	  {"EphemDir", 		required_argument, 0, 	'E'}, 
	  {"DataDir", 		required_argument, 0, 	'D'},
	  {"IFO", 		required_argument, 0, 	'I'},
	  {"SignalOnly", 	no_argument, 0, 	'S'},
	  {"Spin", 		required_argument, 0, 	's'},
	  {"dSpin", 		required_argument, 0, 	'e'},
	  {"SpinBand", 		required_argument, 0, 	'm'},
	  {"EstimSigParam", 	no_argument, 0, 	'p'},
	  {"BaseName",		required_argument,0, 	'i'},
	  {"Fthreshold", 	required_argument, 0, 	'F'},
	  {"useMetric", 	required_argument, 0, 	'M'},
	  {"metricMismatch",	required_argument, 0, 	'X'},
	  {"flipTiling", 	no_argument, 0, 	'x'},
	  {"help", 		no_argument, 0, 	'h'},
	  {"debug", 		required_argument, 0, 	'v'},
	  {0, 0, 0, 0}
	};

      c = getopt_long(argc, argv, "a:b:c:D:d:E:e:F:f:g:hI:i:l:M:m:pr:Ss:t:v:X:xy:z:", long_options, &option_index);

      if (c == -1) 
	break;
      
      switch (c)
	{
	case 'f':
	  /* first search frequency */
	  CLA->Freq=atof(optarg);
	  break;
	case 'b':
	  /* frequency band*/
	  CLA->FreqBand=atof(optarg);
	  break;
	case 'r':
	  /* frequency resolution */
	  CLA->dFreq=atof(optarg);
	  break;
	case 'd':
	  /* sky position delta -- latitude */
	  CLA->Delta=atof(optarg);
	  break;
	case 'g':
	  /* sky position delta resolution -- latitude */
	  CLA->dDelta=atof(optarg);
	  break;
	case 'c':
	  /* sky position delta band -- latitude */
	  CLA->DeltaBand=atof(optarg);
	  break;
	case 'a':
	  /* sky position alpha -- longitude */
	  CLA->Alpha=atof(optarg);
	  break;
	case 'l':
	  /* sky position alpha -- longitude */
	  CLA->dAlpha=atof(optarg);
	  break;
	case 'z':
	  /* sky position alpha -- longitude */
	  CLA->AlphaBand=atof(optarg);
	  break;
	case 't':
	  /* frequency bandwidth */
	  CLA->Dterms=atof(optarg);
	  break;
	case 'y':
	  /* observation year for e-files */
	  CLA->EphemYear=optarg;
	  break;
	case 'E':
	  /* e-file directory */
	  CLA->EphemDir=optarg;
	  break;
	case 'D':
	  /* SFT directory */
	  CLA->DataDir=optarg;
	  break;
	case 'I':
	  /* IFO number */
	  CLA->IFO=atof(optarg);
	  break;
	case 'S':
	  /* frequency bandwidth */
	  CLA->noise=1;
	  break;
	case 's':
	  /* Spin down order */
	  CLA->Spin=atof(optarg);
	  break;
	case 'm':
	  /* Spin down order */
	  CLA->SpinBand=atof(optarg);
	  break;
	case 'p':
	  /* Turn on estimating signal parameters via ML */
	  CLA->EstimSigParam=1;
	  break;
	case 'F':
	  /* Set the threshold for selection of F lines */
	  CLA->Fthreshold=atof(optarg);
	  break;
	case 'e':
	  /* Spin down order */
	  CLA->dSpin=atof(optarg);
	  break;
	case 'M':
	  /* use-metric */
	  GV.useMetric = atoi (optarg);
	  break;

	case 'X':
	  GV.metricMismatch = atof (optarg);
	  break;
	  
	case 'x':
	  GV.flipTiling = 1;
	  break;
	  
	case 'v':
	  lalDebugLevel = atoi (optarg);
	  break;

         case 'i':
        /* Input Base Name */
         CLA->BaseName=optarg;
         break;
	case 'h':
	  /* print usage/help message */
	  fprintf(stdout,"Arguments are (short alternative arguments in brackets):\n");
	  fprintf(stdout,"\t --Freq(-f)\t\tFLOAT\tStarting search frequency in Hz (not set by default)\n");
	  fprintf(stdout,"\t --FreqBand(-b)\t\tFLOAT\tDemodulation frequency band in Hz (set=0.0 by default)\n");
	  fprintf(stdout,"\t --dFreq(-r)\t\tFLOAT\tDemodulation frequency resolution in Hz (set to 1/(8*Tsft*Nsft) by default)\n");
	  fprintf(stdout,"\t --Dterms(-t)\t\tINTEGER\tNumber of terms to keep in Dirichlet kernel sum (default 16)\n");
	  fprintf(stdout,"\t --Alpha(-a)\t\tFLOAT\tSky position alpha (equatorial coordinates) in radians (default 0.0)\n");
	  fprintf(stdout,"\t --AlphaBand(-z)\tFLOAT\tBand in alpha (equatorial coordinates) in radians (default 0.0)\n");
	  fprintf(stdout,"\t --dAlpha(-l)\t\tFLOAT\tResolution in alpha (equatorial coordinates) in radians (default 0.001)\n");
	  fprintf(stdout,"\t --Delta(-d)\t\tFLOAT\tSky position delta (equatorial coordinates) in radians (default 0.0)\n");
	  fprintf(stdout,"\t --DeltaBand(-c)\tFLOAT\tBand in delta (equatorial coordinates) in radians (default 0.0)\n");
	  fprintf(stdout,"\t --dDelta(-g)\t\tFLOAT\tResolution in delta (equatorial coordinates) in radians (default 0.001)\n");
	  fprintf(stdout,"\t --DataDir(-D)\t\tSTRING\tDirectory where SFT's are located (not set by default) \n");
	  fprintf(stdout,"\t --EphemDir(-E)\t\tSTRING\tDirectory where Ephemeris files are located (not set by default) \n");
	  fprintf(stdout,"\t --EphemYear(-y)\tSTRING\tYear of ephemeris files to be used (not set by default) \n");
	  fprintf(stdout,"\t --IFO(-I)\t\tSTRING\tDetector; must be set to 0=GEO, 1=LLO, 2=LHO or 3=Roman Bar (not set by default) \n");
	  fprintf(stdout,"\t --SignalOnly(-S)\t\tSignal only flag; (Default is signal+noise)\n");
	  fprintf(stdout,"\t --Spin(-s)\t\tFLOAT\tStarting spindown parameter (default 0.0) \n");
	  fprintf(stdout,"\t --SpinBand(-m)\t\tFLOAT\tSpindown band (default 0.0)\n");
	  fprintf(stdout,"\t --dSpin(-e)\t\tFLOAT\tSpindown resolution (default 1/(2*Tobs*Tsft*Nsft)\n");
	  fprintf(stdout,"\t --EstimSigParam(-p)\t\tdo Signal Parameter Estimation (Default: no signal parameter estimation)\n");
	  fprintf(stdout,"\t --Fthreshold(-F)\tFLOAT\tSignal Set the threshold for selection of 2F (default 10.0\n");
	  fprintf(stdout,"\t --BaseName(-i)\t\tSTRING\tThe base name of the input  file you want to read.(Default is *SFT* )\n");
	  fprintf(stdout,"\t --useMetric(-M)\tINT2\tUse a metric template grid, with metric type 1 = PtoleMetric, 2 = CoherentMetric\n");
	  fprintf(stdout,"\t --metricMismatch(-X)\tFLOAT\tMaximal mismatch for metric tiling (Default: 0.02)\n");
	  fprintf(stdout,"\t --flipTiling(-x)\t\tUse flipped coordinate order for tiling: {alpha,delta}\n");
	  fprintf(stdout,"\t --debug(-v)\t\tINT2\tSet lalDebugLevel (Default=0)\n");
	  fprintf(stdout,"\t --help(-h)\t\t\tPrint this message.\n");
	  exit(0);
	  break;
	default:
	  /* unrecognized option */
	  errflg++;
	  fprintf(stderr,"Unrecognized option argument %c. Try ComputeFStatistic -h\n",c);
	  exit(1);
      break;
	}
    }

  if(CLA->IFO == -1)
    {
      fprintf(stderr,"No IFO specified; input with -I option.\n");
      fprintf(stderr,"Try ./ComputeFStatistic -h \n");
      return 1;
    }      
  if(CLA->DataDir == "")
    {
      fprintf(stderr,"No SFT directory specified; input directory with -D option.\n");
      fprintf(stderr,"Try ./ComputeFStatistic -h \n");
      return 1;
    }      
  if(CLA->EphemDir == "")
    {
      fprintf(stderr,"No ephemeris data (earth??.dat, sun??.dat) directory specified; input directory with -E option.\n");
      fprintf(stderr,"Try ./ComputeFStatistic -h \n");
      return 1;
    }      
  if(CLA->EphemYear == "")
    {
      fprintf(stderr,"No ephemeris year (earth??.dat, sun??.dat) directory specified; input year with -y option.\n");
      fprintf(stderr,"Try ./ComputeFStatistic -h \n");
      return 1;
    }      
  if(CLA->Freq == 0.0)
    {
      fprintf(stderr,"No search frequency specified; set with -f option.\n");
      fprintf(stderr,"Try ./ComputeFStatistic -h \n");
     return 1;
    }      
  if(CLA->dDelta == 0.0)
    {
      fprintf(stderr,"Value of Delta resolution ( = 0.0) not allowed.\n");
      return 1;
    }      
  if(CLA->dAlpha == 0.0)
    {
      fprintf(stderr,"Value of Alpha resolution ( = 0.0) not allowed.\n");
      return 1;
    }      

  if(CLA->Alpha > LAL_TWOPI)
    {
      fprintf(stderr,"Value of Alpha > 2pi.\n");
      return 1;
    }



  return errflg;
}



/*******************************************************************************/

int CreateDetector(LALDetector *Detector){

/*   LALDetector Detector;  */
  LALFrDetector detector_params;
  LALDetectorType bar;
  LALDetector Detector1;

/*   detector_params=(LALFrDetector )LALMalloc(sizeof(LALFrDetector)); */
 
  bar=LALDETECTORTYPE_CYLBAR;
  strcpy(detector_params.name,"NAUTILUS");
  detector_params.vertexLongitudeRadians=12.67*LAL_PI/180.0;
  detector_params.vertexLatitudeRadians=41.82*LAL_PI/180.0;
  detector_params.vertexElevation=300.0;
  detector_params.xArmAltitudeRadians=0.0;
  detector_params.xArmAzimuthRadians=44.0*LAL_PI/180.0;

  LALCreateDetector(&status,&Detector1,&detector_params,bar);

  *Detector=Detector1;

  return 0;
}

/*******************************************************************************/

int Freemem(void) 
{

  INT4 k;

  /*Free SFTData*/

  for (k=0;k<GV.SFTno;k++)
    {
      LALFree(SFTData[k]->fft->data->data);
      LALFree(SFTData[k]->fft->data);
      LALFree(SFTData[k]->fft);
      LALFree(SFTData[k]);
    }
  LALFree(SFTData);

  /*Free timestamps*/
  LALFree(timestamps);

  LALFree(Fstat.F);
  if(GV.EstimSigParam) 
    {
      LALFree(Fstat.Fa);
      LALFree(Fstat.Fb);
    }

  /*Free DemodParams*/
  LALSDestroyVector(&status, &(DemodParams->amcoe->a));
  LALSDestroyVector(&status, &(DemodParams->amcoe->b));

  LALFree(DemodParams->skyConst);
  LALFree(DemodParams->spinDwn);
  LALFree(DemodParams);

  /* Free DopplerScan-stuff (grid) */
  FreeDopplerScan (&status, &thisScan);

  if (highFLines->clusters) LALFree(highFLines->clusters);
  if (highFLines->Iclust) LALFree(highFLines->Iclust);
  if (highFLines->NclustPoints) LALFree(highFLines->NclustPoints);
    

  /*
  if (highFLines != NULL){
    free(highFLines->NclustPoints);
    free(highFLines->Iclust);
    free(highFLines->clusters);
    free(highFLines);
  }
  */ /* Did not allocate memory for highFLines using malloc/calloc */ 
  LALCheckMemoryLeaks();
  
  return 0;

} /* Freemem() */

/*******************************************************************************/

/* Sorting function to sort into DECREASING order */
int compare(const void *ip, const void *jp)
{
  REAL8 di, dj;

  di=Fstat.F[*(int *)ip];
  dj=Fstat.F[*(int *)jp];

  if (di<dj)
    return 1;
  
  if (di==dj)
    return 0;

  return -1;
}

/*******************************************************************************/

/* This routine prints the values of (f,FF) above a certain threshold */
/* in 2F, called 2Fthr.  If there are more than ReturnMaxN of these, */
/* then it simply returns the top ReturnMaxN of them. If there are */
/* none, then it returns none.  It also returns some basic statisical */
/* information about the distribution of 2F: the mean and standard */
/* deviation. */
/* Returns zero if all is well, else nonzero if a problem was encountered. */
/*    Basic strategy: sort the array by values of F, then look at the */
/*    top ones. Then search for the points above threshold. */


INT4 PrintTopValues(REAL8 TwoFthr, INT4 ReturnMaxN)
{

  INT4 *indexes,i,j,iF,ntop,err,N;
  REAL8 mean=0.0, std=0.0 ,log2 /*=log(2.0)*/;

  log2=medianbias;



  for (i=0;i<highFLines->Nclusters;i++){
    N=highFLines->NclustPoints[i];
    for (j=0;j<N;j++){
      iF=highFLines->Iclust[j];
      Fstat.F[iF]=0.0;
    }/*  end j loop over points of i-th cluster  */
  }/*  end i loop over different clusters */




/*    check that ReturnMaxN is sensible */
  if (ReturnMaxN>GV.FreqImax){
    fprintf(stderr,"PrintTopValues() WARNING: resetting ReturnMaxN=%d to %d\n",
	    ReturnMaxN, GV.FreqImax);
    ReturnMaxN=GV.FreqImax;
  }

/*    create an array of indexes */
  if (!(indexes=(INT4 *)LALMalloc(sizeof(INT4)*GV.FreqImax))){
    fprintf(stderr,"Unable to allocate index array in PrintTopValues()\n");
    return 1;
  }

/*    populate it */
  for (i=0;i<GV.FreqImax;i++)
    indexes[i]=i;

/*   sort array of indexes */
  qsort((void *)indexes, (size_t)GV.FreqImax, sizeof(int), compare);


/*    Normalize */
  TwoFthr*=0.5/log2;

#ifdef FILE_FMAX  
#if 1
/*    print out the top ones */
  for (ntop=0; ntop<ReturnMaxN; ntop++)
    if (Fstat.F[indexes[ntop]]>TwoFthr){
      err=fprintf(fpmax, "%20.10f %10.8f %10.8f %20.15f\n",
		  GV.Freq+indexes[ntop]*GV.dFreq,
		  Alpha, Delta, 2.0*log2*Fstat.F[indexes[ntop]]);
      if (err<=0) {
	fprintf(stderr,"PrintTopValues couldn't print to Fmax!\n");
	LALFree(indexes);
	return 3;
      }
    }
    else
/*        Since array sorted, as soon as value too small we can break */
      break;
#endif
#endif

  /*  find out how many points have been set to zero (N) */
  N=0;
  if (highFLines) {
    for (i=0;i<highFLines->Nclusters; i++){
      N=N+highFLines->NclustPoints[i];
    }
  }
/*    get mean of F[] */
  for (i=0;i<GV.FreqImax; i++)
    mean+=Fstat.F[i];
  mean/=(GV.FreqImax-N);

/*    get sigma for F[] */
  for (i=0; i<GV.FreqImax; i++){
    REAL8 diff=Fstat.F[i]-mean;
    std+=diff*diff;
  }
  std/=(GV.FreqImax-1-N);

/*    normalize by appropriate factors */
  mean*=(2.0*log2);
  std=2.0*log2*sqrt(std);


/*    Find number of values above threshold (could be made O(log N) */
/*    with binary search! */
  for (i=0; i<GV.FreqImax; i++)
    if (Fstat.F[indexes[i]]<=TwoFthr)
      break;

#ifdef FILE_FMAX_DEBG    
/*    print the output */
  err=fprintf(fpmax,"%10.5f %10.8f %10.8f    %d %10.5f %10.5f %10.5f\n",GV.Freq,
	      Alpha, Delta, GV.FreqImax-N, mean, std, 2.0*log2*Fstat.F[indexes[0]]);
#endif
  LALFree(indexes);
#ifdef FILE_FMAX_DEBG    
  if (err<=0) {
    fprintf(stderr,"PrintTopValues couldn't print to Fmax!\n");
    return 4;
  }
#endif

  return 0;
}

/*******************************************************************************/

INT4 EstimatePSDLines(void)
{
#ifdef FILE_PSD
  FILE *outfile;
#endif
#ifdef FILE_PSDLINES
  FILE *outfile1;
#endif
  INT4 i,j,Ntot;                         /* loop indices */
  INT4 nbins=GV.ifmax-GV.ifmin+1;   /* Number of points in SFT's */
  REAL8Vector *Sp=NULL; 
  REAL8Vector *FloorSp=NULL;                        /* Square of SFT */
  INT2 windowSize=100;                  /* Running Median Window Size*/
  REAL4 THR=10000.0;
  
  REAL4 xre,xim,freq;
  REAL8 r0,r1,r2;

  OutliersInput  *outliersInput;
  OutliersParams *outliersParams;
  Outliers       *outliers;
  ClustersInput  *clustersInput;
  ClustersParams *SpClParams;
  Clusters       *SpLines=highSpLines;
    
  INT2 smallBlock=3;
  INT4 wings;

  nbins=(UINT4)nbins;
  wings=windowSize/2;

#ifdef FILE_PSD
  /*  file contains freq, PSD, noise floor */
  if(!(outfile=fopen("PSD.txt","w"))){
    printf("Cannot open PSD.txt file");
    return 1;
  } 
#endif 
#ifdef FILE_PSDLINES
  /*  file contains freq, PSD, noise floor,lines */
  if(!(outfile1=fopen("PSDLines.txt","w"))){
    printf("Cannot open PSD.txt file");
    return 1;
  }
#endif

  /* Allocate memory for input & output */
  /* if (!(Sp = (double *) calloc(nbins,sizeof(double)))){ */
  /*   printf("Memory allocation failure"); */
  /*   return 0; */
  /* } */
  
  LALDCreateVector(&status, &Sp, nbins);
  LALDCreateVector(&status, &FloorSp, nbins);
  
  
  /* loop over each SFTs */
  for (i=0;i<GV.SFTno;i++){
    
    /* loop over SFT data to estimate noise */
    for (j=0;j<nbins;j++){
      xre=SFTData[i]->fft->data->data[j].re;
      xim=SFTData[i]->fft->data->data[j].im;
      Sp->data[j]=Sp->data[j]+(REAL8)(xre*xre+xim*xim);
    }
  }/*end loop over SFTs*/
  
  /*Average Sp*/
  for (j=0;j<nbins;j++){
    Sp->data[j]=Sp->data[j]/GV.SFTno;
  }
  Sp->length=nbins;
  FloorSp->length=nbins;

  j=EstimateFloor(Sp, windowSize, FloorSp);
 
  if (!(outliers=(Outliers *)LALMalloc(sizeof(Outliers)))){
    fprintf(stderr,"Memory allocation failure for SpOutliers\n");
    return 1;
  }
  outliers->Noutliers=0;

  if (!(outliersParams=(OutliersParams *)LALMalloc(sizeof(OutliersParams)))){
    fprintf(stderr,"Memory allocation failure for OutliersParams\n");
    return 1;
  }
  if (!(outliersInput=(OutliersInput *)LALMalloc(sizeof(OutliersInput)))){
    fprintf(stderr,"Memory allocation failure for OutliersParams\n");
    return 1;
  }
  
  outliersParams->Thr=THR;
  outliersParams->Floor = FloorSp;
  outliersParams->wings=wings; /*these must be the same as ClustersParams->wings */
  outliersInput->ifmin=GV.ifmin;
  outliersInput->data = Sp;

  ComputeOutliers(outliersInput, outliersParams, outliers);

   if (outliers->Noutliers == 0){

#ifdef FILE_PSD
     /*  PSD.txt file contains freq, PSD, noise floor   */
     for (i=0;i<nbins;i++){ 
       freq=(GV.ifmin+i)/GV.tsft;
       r0=Sp->data[i];
       r1=FloorSp->data[i];
       fprintf(outfile,"%f %E %E\n",freq,r0,r1);
     }
#endif

#ifdef FILE_PSD     
     fclose(outfile);
#endif
#ifdef FILE_PSDLINES
     fclose(outfile1);
#endif

     LALFree(outliers->ratio);
     LALFree(outliers);
     LALFree(outliersParams);
     LALFree(outliersInput);
     LALDDestroyVector(&status,&Sp);
     LALDDestroyVector(&status,&FloorSp);

     return 0;

   }
  


   if (!(SpClParams=(ClustersParams *)LALMalloc(sizeof(ClustersParams)))){ 
     printf("Memory allocation failure for SpClusterParams");
     return 1;
   }

   if (!(clustersInput=(ClustersInput *)LALMalloc(sizeof(ClustersInput)))){ 
     printf("Memory allocation failure for SpClusters");
     return 1;
   }
      
   SpClParams->wings=wings;
   SpClParams->smallBlock=smallBlock;
   
   clustersInput->outliersInput = outliersInput;
   clustersInput->outliersParams= outliersParams;
   clustersInput->outliers      = outliers;     
   
   j=DetectClusters(clustersInput, SpClParams, SpLines);
   if (j!=0){
     printf("DetectClusters problem");
     return 1;
   }
      
   /*  sum of points in all lines */
   Ntot=0;
   for (i=0;i<SpLines->Nclusters;i++){ 
     Ntot=Ntot+SpLines->NclustPoints[i];
   }
   
#ifdef FILE_PSDLINES
   /*  PSDLines file contains: PSD, noise floor and lines. */
   for (i=0;i<Ntot;i++){ 
     j=SpLines->Iclust[i];
     freq=(GV.ifmin+SpLines->Iclust[i])/GV.tsft;
     r0=Sp->data[j];
     r1=FloorSp->data[j];
     r2=SpLines->clusters[i]*FloorSp->data[j];
     fprintf(outfile1,"%f %E %E %E\n",freq,r0,r1,r2);
   }
#endif

#ifdef FILE_PSD   
   /*  PSD.txt file contains freq, PSD, noise floor   */
   for (i=0;i<nbins;i++){ 
     freq=(GV.ifmin+i)/GV.tsft;
     r0=Sp->data[i];
     r1=FloorSp->data[i];
     fprintf(outfile,"%f %E %E\n",freq,r0,r1);
   }
#endif

#ifdef FILE_PSD   
   fclose(outfile);
#endif
#ifdef FILE_PSDLINES
   fclose(outfile1);
#endif

   LALFree(outliers->ratio);
   LALFree(outliers->outlierIndexes);
   LALFree(outliers);
   LALFree(outliersParams);
   LALFree(outliersInput);
   LALDDestroyVector(&status,&Sp);
   LALDDestroyVector(&status,&FloorSp);
   LALFree(SpClParams);
   LALFree(clustersInput);

   return 0;
}

/*******************************************************************************/

INT4 EstimateFLines(void)
{
#ifdef FILE_FTXT  
  FILE *outfile;
#endif
#ifdef FILE_FLINES  
  FILE *outfile1;
#endif
  INT4 i,j,Ntot;                         /* loop indices */
  INT4 nbins=GV.FreqImax;                /* Number of points in F */
  REAL8Vector *F1=NULL; 
  REAL8Vector *FloorF1=NULL;                        /* Square of SFT */
  /* INT2 windowSize=(0.01/GV.dFreq);               0.1 is 1E-4*1000 */
  INT2 windowSize=100;
  REAL4 THR=10.0;
  
  REAL4 freq;
  REAL8 r0,r1,r2;

  REAL8 dmp;

  OutliersInput  *outliersInput;
  OutliersParams *outliersParams;
  Outliers       *outliers;
  ClustersInput  *clustersInput;
  ClustersParams *SpClParams;
  Clusters       *SpLines=highFLines;
    
  INT2 smallBlock=1;
  
  INT4 wings;

  nbins=(UINT4)nbins;

  THR=GV.Fthreshold;


  /* wings=windowSize/2; */
  /*  0.0002 is the max expected width of the F stat curve for signal */
  /*  with ~ 10 h observation time */
  /*  0.0001 = 0.0002/2 */
  /*  let me put 0.005 */
  dmp=0.5+0.0002/GV.dFreq;
  wings=dmp;


  if (windowSize > nbins){
    windowSize = nbins/2.0;
    /* printf("Had to change windowSize for running median in F floor estimate\n"); */
  }

#ifdef FILE_FTXT
  /*  file contains freq, PSD, noise floor */
  if(!(outfile=fopen("F.txt","w"))){
    printf("Cannot open F.txt file\n");
    return 1;
  }
#endif
#ifdef FILE_FLINES  
  /*  file contains freq, PSD, noise floor,lines */
  if(!(outfile1=fopen("FLines.txt","w"))){
    printf("Cannot open FLines.txt file\n");
    return 1;
  }
#endif


  LALDCreateVector(&status, &F1, nbins);
  LALDCreateVector(&status, &FloorF1, nbins);
    
  /* loop over SFT data to estimate noise */
  for (j=0;j<nbins;j++){
    F1->data[j]=Fstat.F[j];
    FloorF1->data[j]=1.0;
  }
  
  F1->length=nbins;
  FloorF1->length=nbins;

  /*   j=EstimateFloor(F1, windowSize, FloorF1); */
 
  if (!(outliers=(Outliers *)LALMalloc(sizeof(Outliers)))){
    fprintf(stderr,"Memory allocation failure for SpOutliers\n");
    return 1;
  }
  outliers->Noutliers=0;

  if (!(outliersParams=(OutliersParams *)LALMalloc(sizeof(OutliersParams)))){
    fprintf(stderr,"Memory allocation failure for OutliersParams\n");
    return 1;
  }
  if (!(outliersInput=(OutliersInput *)LALMalloc(sizeof(OutliersInput)))){
    fprintf(stderr,"Memory allocation failure for OutliersParams\n");
    return 1;
  }
  
  outliersParams->Thr=THR/(2.0*medianbias);
  outliersParams->Floor = FloorF1;
  outliersParams->wings=wings; /*these must be the same as ClustersParams->wings */
  outliersInput->ifmin=((GV.Freq/GV.dFreq)+0.5);
  outliersInput->data = F1;

  ComputeOutliers(outliersInput, outliersParams, outliers);

   if (outliers->Noutliers == 0){

#ifdef FILE_FTXT
     /*  F.txt file contains freq, F, noise floor of F   */
     for (i=0;i<nbins;i++){ 
       freq=GV.Freq+i*GV.dFreq;
       r0=F1->data[i];
       r1=FloorF1->data[i];
       fprintf(outfile,"%f %E %E\n",freq,r0,r1);
     }
#endif     

#ifdef FILE_FTXT
     fclose(outfile);
#endif
#ifdef FILE_FLINES  
     fclose(outfile1);
#endif 
     LALFree(outliers->ratio);
     LALFree(outliers);
     LALFree(outliersParams);
     LALFree(outliersInput);
     LALDDestroyVector(&status,&F1);
     LALDDestroyVector(&status,&FloorF1);

     /*      fprintf(stderr,"Nclusters zero \n"); */
     /*      fflush(stderr); */

     return 0;

   }
  
   /*      fprintf(stderr,"Nclusters non zero \n"); */
   /*      fflush(stderr); */

   if (!(SpClParams=(ClustersParams *)LALMalloc(sizeof(ClustersParams)))){ 
     printf("Memory allocation failure for SpClusterParams");
     return 1;
   }

   
   if (!(clustersInput=(ClustersInput *)LALMalloc(sizeof(ClustersInput)))){ 
     printf("Memory allocation failure for SpClusters");
     return 1;
   }
      
   SpClParams->wings=wings;
   SpClParams->smallBlock=smallBlock;
   
   clustersInput->outliersInput = outliersInput;
   clustersInput->outliersParams= outliersParams;
   clustersInput->outliers      = outliers;     
   
   j=DetectClusters(clustersInput, SpClParams, SpLines);
   if (j!=0){
     printf("DetectClusters problem");
     return 1;
   }
   
   
   /*  sum of points in all lines */
   Ntot=0;
   for (i=0;i<SpLines->Nclusters;i++){ 
     Ntot=Ntot+SpLines->NclustPoints[i];
   }
   


#ifdef FILE_FLINES  
   /*  FLines file contains: F, noise floor and lines. */
   for (i=0;i<Ntot;i++){ 
     j=SpLines->Iclust[i];
     freq=(GV.Freq+SpLines->Iclust[i]*GV.dFreq);
     r0=F1->data[j];
     r1=FloorF1->data[j];
     r2=SpLines->clusters[i]*FloorF1->data[j];
     fprintf(outfile1,"%f %E %E %E\n",freq,r0,r1,r2);
   }
#endif
#ifdef FILE_FTXT   
   /*  PSD.txt file contains freq, PSD, noise floor   */
   for (i=0;i<nbins;i++){ 
     freq=GV.Freq+i*GV.dFreq;
     r0=F1->data[i];
     r1=FloorF1->data[i];
     fprintf(outfile,"%f %E %E\n",freq,r0,r1);
   }
#endif   

#ifdef FILE_FTXT
   fclose(outfile);
#endif
#ifdef FILE_FLINES  
   fclose(outfile1);
#endif   

   LALFree(outliers->ratio);
   LALFree(outliers->outlierIndexes);
   LALFree(outliers);
   LALFree(outliersParams);
   LALFree(outliersInput);
   LALDDestroyVector(&status,&F1);
   LALDDestroyVector(&status,&FloorF1);
   LALFree(SpClParams);
   LALFree(clustersInput);

   return 0;


}

/*******************************************************************************/
/*******************************************************************************/

INT4 NormaliseSFTDataRngMdn(void)
{
#ifdef FILE_SPRNG  
  FILE *outfile;
#endif
  INT4 i,j,m,lpc,il;                         /* loop indices */
  INT4 Ntot,nbins=GV.ifmax-GV.ifmin+1;   /* Number of points in SFT's */
  REAL8Vector *Sp=NULL, *RngMdnSp=NULL;   /* |SFT|^2 and its rngmdn  */
  REAL8 B;                          /* SFT Bandwidth */
  REAL8 deltaT,norm,*N, *Sp1;
  INT2 windowSize=50;                  /* Running Median Window Size*/
  REAL4 xre,xim,xreNorm,ximNorm;


  /* The running median windowSize in this routine determines 
     the sample bias which, instead of log(2.0), must be 
     multiplied by F statistics.
  */

  if (GV.noise != 1)
    MedianBias(&windowSize,&medianbias);

  LALDCreateVector(&status, &Sp, (UINT4)nbins);
  LALDCreateVector(&status, &RngMdnSp, (UINT4)nbins);

  nbins=(INT2)nbins;

  if(!(N= (REAL8 *) LALCalloc(nbins,sizeof(REAL8)))){ 
    printf("Memory allocation failure");
    return 0;
  }
   if(!(Sp1= (REAL8 *) LALCalloc(nbins,sizeof(REAL8)))){ 
    printf("Memory allocation failure");
    return 0;
  }
   
  
  /* loop over each SFTs */
  for (i=0;i<GV.SFTno;i++)         
    {
      /* Set to zero the values */
      for (j=0;j<nbins;j++){
	RngMdnSp->data[j] = 0.0;
	Sp->data[j]       = 0.0;
      }
      
      /* loop over SFT data to estimate noise */
      for (j=0;j<nbins;j++){
	xre=SFTData[i]->fft->data->data[j].re;
	xim=SFTData[i]->fft->data->data[j].im;
	Sp->data[j]=(REAL8)(xre*xre+xim*xim);
      }
      
      /* Compute running median */
      EstimateFloor(Sp, windowSize, RngMdnSp);
      
      /* compute how many cluster points in all */
      /* substitute the line profiles value in RngMdnSp */
      Ntot=0;
      if (highSpLines != NULL){
	for (il=0;il<highSpLines->Nclusters;il++){
	  Ntot=Ntot+highSpLines->NclustPoints[il];
	}
	for (j=0;j<Ntot;j++){
	  m=highSpLines->Iclust[j];
	  RngMdnSp->data[m]=RngMdnSp->data[m]*highSpLines->clusters[j];	
	}
      }

     
      /*Compute Normalization factor*/
      /* for signal only case as well */  
      for (lpc=0;lpc<nbins;lpc++){
	N[lpc]=1.0/sqrt(2.0*RngMdnSp->data[lpc]);
      }
      
      if(GV.noise == 1){
	B=(1.0*GV.nsamples)/(1.0*GV.tsft);
	deltaT=1.0/(2.0*B);
	norm=deltaT/sqrt(GV.tsft);
	for (lpc=0;lpc<nbins;lpc++){
	  N[lpc]=norm;
	}
      }
      
      /*  loop over SFT data to normalise it (with N) */
      /*  also compute Sp1, average normalized PSD */
      /*  and the sum of the PSD in the band, SpSum */
      for (j=0;j<nbins;j++){
	xre=SFTData[i]->fft->data->data[j].re;
	xim=SFTData[i]->fft->data->data[j].im;
	xreNorm=N[j]*xre; 
	ximNorm=N[j]*xim; 
	SFTData[i]->fft->data->data[j].re = xreNorm;    
	SFTData[i]->fft->data->data[j].im = ximNorm;
	Sp1[j]=Sp1[j]+xreNorm*xreNorm+ximNorm*ximNorm;
      }
      
    } /* end loop over SFTs*/

#ifdef FILE_SPRNG  
  if(!(outfile=fopen("SpRng.txt","w"))){ 
    printf("Cannot open output file"); 
    return 1;
  } 


  for (j=0;j<nbins;j++){
    Sp1[j]=2.0*Sp1[j]/(1.0*GV.SFTno);
    fprintf(outfile,"%f %E \n",(GV.ifmin+j)/GV.tsft,Sp1[j]); 
  }
  
  fclose(outfile);
#endif

  
  LALFree(N);
  LALFree(Sp1);
  LALDDestroyVector(&status, &RngMdnSp);
  LALDDestroyVector(&status, &Sp);
  
  return 0;
}

/*******************************************************************************/
