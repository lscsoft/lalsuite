/*********************************************************************************/
/*                    F-statistic generation code for known pulsars              */
/*                                                                               */
/*			      Y. Ioth, M.A. Papa, X. Siemens, R.Prix             */
/*                                                                               */
/*                 Albert Einstein Institute/UWM - started September 2002        */
/*********************************************************************************/
#include <lal/LALDemod.h>
#include <lal/RngMedBias.h>

#include "ComputeFStatistic.h"
#include "rngmed.h"
#include "clusters.h"
#include "DopplerScan.h"

#include "UserInput.h"

NRCSID( COMPUTEFSTATISTIC, "$Id$");

/*
#define DEBG_FAFB                
#define DEBG_ESTSIGPAR
#define DEBG_MAIN 
#define DEBG_SGV 
*/

/* If FILE_FILENAME is defined, then print the corresponding file  */
#define FILE_FSTATS
#define FILE_FMAX

/*
#define FILE_FLINES
#define FILE_FTXT 
#define FILE_PSD 
#define FILE_PSDLINES 
#define FILE_SPRNG 
*/

/*----------------------------------------------------------------------*/
/* 
 * User-specified variables: provided either by default, config-file or command-line
 * NOTE: these HAVE to follow the name-convention "uvar_VarName"
 */

INT4 uvar_Dterms 	= 16;
INT4 uvar_IFO 		= -1;
BOOLEAN uvar_SignalOnly = 0;
BOOLEAN uvar_EstimSigParam = 0;
REAL8 uvar_Freq 	= 0.0;
REAL8 uvar_dFreq 	= 0.0;
REAL8 uvar_FreqBand 	= 0.0;
REAL8 uvar_Alpha 	= 0.0;
REAL8 uvar_dAlpha 	= 0.001;
REAL8 uvar_AlphaBand 	= 0.0;
REAL8 uvar_Delta 	= 0.0;
REAL8 uvar_dDelta 	= 0.001;
REAL8 uvar_DeltaBand 	= 0.0;
REAL8 uvar_Spin;
REAL8 uvar_dSpin 	= 0.0;
REAL8 uvar_SpinBand;
REAL8 uvar_Fthreshold 	= 10.0;
INT4 uvar_EphemYear	= -1;
INT4 uvar_Metric 	=  LAL_METRIC_NONE;	
REAL8 uvar_metricMismatch = 0.02;
CHAR *uvar_skyRegion	= NULL;
CHAR *uvar_DataDir	= NULL;
CHAR *uvar_EphemDir	= NULL;
CHAR *uvar_BaseName	= NULL;
BOOLEAN uvar_debug;		/* NOTE: don't use this, it's a dummy!! lalDebugLevel is set independently */
BOOLEAN uvar_help;
CHAR *uvar_outputLabel 	= NULL;

UserVariable uvars[] = {
  regUserVar (Dterms, 	UVAR_INT4, 't', "Number of terms to keep in Dirichlet kernel sum"),
  regUserVar (Freq, 	UVAR_REAL8,'f', "Starting search frequency in Hz"),
  regUserVar (FreqBand, UVAR_REAL8,'b', "Demodulation frequency band in Hz"),
  regUserVar (dFreq, 	UVAR_REAL8,'r', "Demodulation frequency resolution in Hz (set to 1/(8*Tsft*Nsft) by default"),
  regUserVar (Alpha, 	UVAR_REAL8,'a', "Sky position alpha (equatorial coordinates) in radians"),
  regUserVar (AlphaBand,UVAR_REAL8,'z', "Band in alpha (equatorial coordinates) in radians"),
  regUserVar (dAlpha, 	UVAR_REAL8,'l', "Resolution in alpha (equatorial coordinates) in radians"),
  regUserVar (Delta, 	UVAR_REAL8,'d', "Sky position delta (equatorial coordinates) in radians"),
  regUserVar (DeltaBand,UVAR_REAL8,'c', "Band in delta (equatorial coordinates) in radians"),
  regUserVar (dDelta, 	UVAR_REAL8,'g', "Resolution in delta (equatorial coordinates) in radians"),
  regUserVar (DataDir, 	UVAR_CHAR, 'D', "Directory where SFT's are located"),
  regUserVar (EphemDir, UVAR_CHAR, 'E', "Directory where Ephemeris files are located"),
  regUserVar (EphemYear,UVAR_INT4, 'y', "Year of ephemeris files to be used"),
  regUserVar (IFO, 	UVAR_INT4, 'I', "Detector, must be set to 0=GEO, 1=LLO, 2=LHO or 3=Roman Bar"),
  regUserVar (SignalOnly,UVAR_BOOL,'S', "Signal only flag"),
  regUserVar (Spin, 	UVAR_REAL8,'s', "Starting spindown parameter"),
  regUserVar (SpinBand, UVAR_REAL8,'m', "Spindown band"),
  regUserVar (dSpin, 	UVAR_REAL8,'e', "Spindown resolution (default 1/(2*Tobs*Tsft*Nsft)"),
  regUserVar (EstimSigParam, UVAR_BOOL, 'p', "Do Signal Parameter Estimation"),
  regUserVar (Fthreshold,UVAR_REAL8,'F', "Signal Set the threshold for selection of 2F"),
  regUserVar (BaseName, UVAR_CHAR, 'i', "The base name of the input  file you want to read"),
  regUserVar (Metric,	UVAR_INT4, 'M', "Metric for template grid: 0=none, 1 = PtoleMetric, 2 = CoherentMetric"),
  regUserVar (metricMismatch, UVAR_REAL8, 'X', "Maximal mismatch for metric tiling"),
  regUserVar (debug, 	UVAR_INT4, 'v', "Set lalDebugLevel"),
  regUserVar (help, 	UVAR_BOOL, 'h', "Print this message"),
  regUserVar (skyRegion,UVAR_CHAR, 'R', "Specify sky-region by polygon"),
  regUserVar (outputLabel,UVAR_CHAR,'o', "Label to be appended to all output file-names"),
  {NULL, 0, 0, NULL, NULL}
};


/*----------------------------------------------------------------------*/


FFT **SFTData=NULL;                 /* SFT Data for LALDemod */
DemodPar *DemodParams  = NULL;      /* Demodulation parameters for LALDemod */
LIGOTimeGPS *timestamps=NULL;       /* Time stamps from SFT data */
LALFstat Fstat;
INT4 lalDebugLevel=0,i,a,d,s,irec;
static LALStatus status;
AMCoeffs amc;
REAL8 MeanOneOverSh=0.0;
REAL8 Alpha,Delta;
Clusters HFLines, HPLines;
Clusters *highSpLines=&HPLines, *highFLines=&HFLines;
/* #ifdef FILE_FMAX     */
FILE *fpmax;
/* #endif */
/* #ifdef FILE_STATS */
FILE *fpstat;
/* #endif     */
REAL8 medianbias=1.0;

DopplerScanState thisScan;
ConfigVariables GV;

/* local prototypes */
INT2 SetGlobalVariables(ConfigVariables *cfg);
INT4 ReadSFTData(void);
INT4 CreateDemodParams(void);
INT4 NormaliseSFTData(void);
INT4 CreateDetector(LALDetector *Detector);
INT4 AllocateMem(void);
INT4 Freemem(void);

INT4 EstimateSignalParameters(INT4 * maxIndex);
INT4 EstimatePSDLines(void);
INT4 EstimateFLines(void);
INT4 writeFLines(INT4 *maxIndex);
INT4 NormaliseSFTDataRngMdn(void);
INT4 PrintTopValues(REAL8 TwoFthr, INT4 ReturnMaxN);
INT4 EstimateFloor(REAL8Vector *Sp, INT2 windowSize, REAL8Vector *SpFloor);
int compare(const void *ip, const void *jp);
INT4 writeFaFb(INT4 *maxIndex);

/*----------------------------------------------------------------------
 * MAIN
 *----------------------------------------------------------------------*/
int main(int argc,char *argv[]) 
{
  INT4 *maxIndex=NULL; /*  array that contains indexes of maximum of each cluster */
  DopplerPosition dopplerpos;
  DopplerScanInit scanInit;
  CHAR Fstatsfilename[256];         /* Fstats file name*/
  CHAR Fmaxfilename[256];           /* Fmax file name*/

  /*----------------------------------------------------------------------
   * ok, this is awkward, but the LALDebugLevel has to be set _first thing_, 
   * before any calls to LALMalloc etc have been made, so we do it my hand: 
   */
  for (i=1; i < argc; i++) {
    if ( strncmp (argv[i], "-v", 2) == 0 ) {
      if (strlen(argv[i]) == 2) {
	if (i+1<argc) lalDebugLevel = atoi (argv[i+1]);
      } else { lalDebugLevel = atoi (argv[i]+2); }
    }
  } /* for cmd-line */
  /*----------------------------------------------------------------------*/

  LALReadUserInput (&status, argc,argv, uvars);

  /* print help-string and exit if -h was specified */
  if (uvar_help)
    {
      /* little hack here: default for uvar_help = False, therefore we have to set that back: */
      CHAR *helpstr = NULL;
      uvar_help = 0;
      LALGetUvarHelpString (&status, &helpstr, uvars);
      printf ("Arguments are (short alternative arguments in brackets):\n");
      printf (helpstr);
      LALFree (helpstr);
      exit (0);
    }

  if (SetGlobalVariables (&GV) ) return 2;

  if (AllocateMem()) return 3;

  if (ReadSFTData()) return 4;

  /*  This fills-in highSpLines that are then used by NormaliseSFTRngMdn */
#if 0
  if (GV.SignalOnly!=1){
    if (EstimatePSDLines()) return 6;
  }
#endif

  if (NormaliseSFTDataRngMdn()) return 5;

#ifdef FILE_FMAX  
  /*   open file */
  strcpy(Fmaxfilename,"Fmax");
  if (uvar_outputLabel)
    strcat(Fmaxfilename,uvar_outputLabel);
  if (!(fpmax=fopen(Fmaxfilename,"w"))){
    fprintf(stderr,"in Main: unable to open Fmax file\n");
    return 2;
  }
#endif
#ifdef FILE_FSTATS  
  /*      open file */
  strcpy(Fstatsfilename,"Fstats");
  if (uvar_outputLabel)
    strcat(Fstatsfilename,uvar_outputLabel);
  if (!(fpstat=fopen(Fstatsfilename,"w"))){
    fprintf(stderr,"in Main: unable to open Fstats file\n");
    return 2;
  }
#endif
	
  scanInit.dAlpha = uvar_dAlpha;
  scanInit.dDelta = uvar_dDelta;

  /* use metric-grid? and if so, for what maximal mismatch? */
  scanInit.metricType = uvar_Metric;
  scanInit.metricMismatch = uvar_metricMismatch;

  scanInit.obsBegin = SFTData[0]->fft->epoch;
  scanInit.obsDuration = SFTData[GV.SFTno-1]->fft->epoch.gpsSeconds - scanInit.obsBegin.gpsSeconds + GV.tsft;

  /*  scanInit.obsDuration = 36000; */

  scanInit.fmax  = uvar_Freq;
  if (uvar_FreqBand > 0) scanInit.fmax += uvar_FreqBand;
  scanInit.Detector = GV.Detector;

  if (uvar_skyRegion)
    {
      scanInit.skyRegion = LALMalloc (strlen (uvar_skyRegion) + 1);
      strcpy (scanInit.skyRegion, uvar_skyRegion);
    }
  else	 /* parse (Alpha+AlphaBand, Delta+DeltaBand) into a sky-region string */
    {
      REAL8 a, d, Da, Dd;
      REAL8 eps = 1.0e-12;	/* slightly push outside the upper and right border (for backwards compatibility) */
      a = uvar_Alpha;
      d = uvar_Delta;
      Da = uvar_AlphaBand + eps;
      Dd = uvar_DeltaBand + eps;

      scanInit.skyRegion = LALMalloc (512); /* should be enough for 4 points... */
      sprintf (scanInit.skyRegion, "(%.16f, %.16f), (%.16f, %.16f), (%.16f, %.16f), (%.16f, %.16f)", 
	      a, d, 
	      a + Da, d, 
	      a + Da, d + Dd,
	      a, d + Dd );
    }

  InitDopplerScan ( &status, &thisScan, scanInit);
  
  LALFree (scanInit.skyRegion);
  
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
      /* Have we scanned all DopplerPositions yet? */
      if (dopplerpos.finished)  
	break;


      Alpha = dopplerpos.skypos.longitude;
      Delta = dopplerpos.skypos.latitude; 
      
      
      if (CreateDemodParams()) return 6;
      /* loop over spin params */
      for(s=0;s<GV.SpinImax;s++)
	{
	  DemodParams->spinDwn[0]=uvar_Spin + s*uvar_dSpin;
	  LALDemod (&status, &Fstat, SFTData, DemodParams);
	  
	  /*  This fills-in highFLines that are then used by discardFLines */
	  if (GV.FreqImax > 5)
	    if (EstimateFLines()) return 6;
	  

#ifdef DEBG_MAIN
	  for(i=0;i < GV.FreqImax ;i++)
	    {
	      /* medianbias is 1 if GV.SignalOnly=1 */ 
	      fprintf(stdout,"%20.10f %e %20.17f %20.17f %20.17f\n",
		      uvar_Freq+i*uvar_dFreq,  DemodParams->spinDwn[0],
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
	  
	  if( uvar_EstimSigParam &&(highFLines !=NULL) && (highFLines->Nclusters >0))
	    if(writeFaFb(maxIndex)) return 255;
	  
	  
	  if( uvar_EstimSigParam &&(highFLines !=NULL) &&(highFLines->Nclusters >0)) {
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
  CHAR Paramfilename[256];
  

#ifdef DEBG_ESTSIGPAR
  REAL8 Ftest;
  REAL8 A2test,A3test,A4test;
#endif


  strcpy(Paramfilename,"ParamMLE");
  if (uvar_outputLabel)
    strcat(Paramfilename,uvar_outputLabel);
  
  if(!(fpMLEParam=fopen(Paramfilename,"w")))
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
		irec,uvar_Freq+irec*uvar_dFreq,A1,A1test);
	fprintf(stderr,"relative error Abs((A1-A1test)/A1)=%lf\n",
		fabs(A1-A1test)/fabs(A1));
	exit(1);
      }
      if(fabs(A2-A2test)>fabs(A2)/(10e5)){ 
	fprintf(stderr,"Something is wrong with Estimate A2\n");
	fprintf(stderr,"Frequency index %d, %lf (Hz),A2=%f,A2test=%f\n",
		irec,uvar_Freq+irec*uvar_dFreq,A2,A2test);
	fprintf(stderr,"relative error Abs((A2-A2test)/A2)=%lf\n",
		fabs(A2-A2test)/fabs(A2));
	exit(1);
      }
      if(fabs(A3-A3test)>fabs(A3)/(10e5)){ 
	fprintf(stderr,"Something is wrong with Estimate A3\n");
	fprintf(stderr,"Frequency index %d, %lf (Hz),A3=%f,A3test=%f\n",
		irec,uvar_Freq+irec*uvar_dFreq,A3,A3test);
	fprintf(stderr,"relative error Abs((A3-A3test)/A3)=%lf\n",
		fabs(A3-A3test)/fabs(A3));
	exit(1);
      }
      if(fabs(A4-A4test)>fabs(A4)/(10e5)){ 
	fprintf(stderr,"Something is wrong with Estimate A4\n");
	fprintf(stderr,"Frequency index %d, %lf (Hz),A4=%f,A4test=%f\n",
		irec,uvar_Freq+irec*uvar_dFreq,A1,A1test);
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
		irec,uvar_Freq+irec*uvar_dFreq,Fstat.F[irec],Ftest);
	fprintf(stderr,"relative error Abs((F-Ftest)/Ftest)=%lf\n",
		fabs(Fstat.F[irec]-Ftest)/fabs(Ftest));
	exit(1);
      }
#endif


      /* normalization */
      h0mle=h0mle*norm;


      /* For the real data, we need to multiply long(2.0) */
      /* Because we use running median to estimate the S_h. */
      /* if(GV.SignalOnly!=1) 
	h0mle=h0mle*sqrt(medianbias);
      */
      /* medianbias is 1 when GV.SignalOnly==1 */
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
	      A,B,C,uvar_Freq+irec*uvar_dFreq,h0mle,Fstat.F[irec]*medianbias);
      }
#endif

      /* Note that we print out MLE of 2.0*Phi0_JKS */
      /* because Phi0_PULGROUPDOC=2.0*Phi0_JKS */
      /* and Phi0_PULGROUPDOC is the one used in In.data. */
 
      /* medianbias is 1 if GV.SignalOnly==1 */
      fprintf(fpMLEParam,"%16.8lf %22E", uvar_Freq + irec*uvar_dFreq, 2.0*medianbias*Fstat.F[irec]);


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
  CHAR FaFbfilename[256];
  
  strcpy(FaFbfilename,"FaFb");
  if (uvar_outputLabel)
    strcat(FaFbfilename,uvar_outputLabel);
  sprintf(noiseswitch,"%02d", uvar_SignalOnly);
  strcat(FaFbfilename,noiseswitch);

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
	    uvar_Freq+maxIndex[irec]*uvar_dFreq,
	    Fstat.F[maxIndex[irec]]*bias*bias);
    fprintf(fp,"%22.12f %22.12f\n",uvar_Freq + index * uvar_dFreq, uvar_dFreq);
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
	      uvar_Freq+index*uvar_dFreq,Fstat.F[index]*bias*bias,
	      DemodParams->spinDwn[0], Alpha, Delta,
	      Fstat.Fa[index].re/sqrt(GV.SFTno)*bias,
	      Fstat.Fa[index].im/sqrt(GV.SFTno)*bias,
	      Fstat.Fb[index].re/sqrt(GV.SFTno)*bias,
	      Fstat.Fb[index].im/sqrt(GV.SFTno)*bias,
	      amc.A,amc.B,amc.C);
#else
      /* Freqency, Re[Fa],Im[Fa],Re[Fb],Im[Fb], F */
      fprintf(fp,"%22.16f %22.12f %22.12f %22.12f %22.12f %22.12f\n",
	      uvar_Freq+index*uvar_dFreq,
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

  DemodParams->f0=uvar_Freq;
  DemodParams->imax=GV.FreqImax;
  DemodParams->df=uvar_dFreq;

  DemodParams->Dterms=uvar_Dterms;
  DemodParams->ifmin=GV.ifmin;

  DemodParams->returnFaFb = uvar_EstimSigParam;

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
    fr=uvar_Freq + imax*uvar_dFreq;
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
      if(uvar_SignalOnly == 1)
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
  INT4 fileno=0,offset;
  FILE *fp;
  size_t errorcode;
  UINT4 ndeltaf;

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

INT2 
SetGlobalVariables(ConfigVariables *cfg)
{

  CHAR command[512];
  CHAR tmp[50];
  FILE *fp;
  size_t errorcode;
  REAL8 df;                         /* freq resolution */
  INT4 fileno=0;   
  glob_t globbuf;

  /* do some sanity checks on the user-input before we proceed */
  if(uvar_IFO == -1)
    {
      fprintf(stderr,"No IFO specified; input with -I option.\n");
      fprintf(stderr,"Try ./ComputeFStatistic -h \n");
      exit(-1);
    }      
  if(uvar_DataDir == NULL)
    {
      fprintf(stderr,"No SFT directory specified; input directory with -D option.\n");
      fprintf(stderr,"Try ./ComputeFStatistic -h \n");
      exit(-1);
    }      
  if(uvar_EphemDir == NULL)
    {
      fprintf(stderr,"No ephemeris data (earth??.dat, sun??.dat) directory specified; input directory with -E option.\n");
      fprintf(stderr,"Try ./ComputeFStatistic -h \n");
      exit(-1);
    }      
  if(uvar_EphemYear == -1)
    {
      fprintf(stderr,"No ephemeris year (earth??.dat, sun??.dat) directory specified; input year with -y option.\n");
      fprintf(stderr,"Try ./ComputeFStatistic -h \n");
      exit(-1);
    }      
  if(uvar_Freq == 0.0)
    {
      fprintf(stderr,"No search frequency specified; set with -f option.\n");
      fprintf(stderr,"Try ./ComputeFStatistic -h \n");
      exit(-1);
    }      
  if(uvar_dDelta == 0.0)
    {
      fprintf(stderr,"Value of Delta resolution ( = 0.0) not allowed.\n");
      exit(-1);
    }      
  if(uvar_dAlpha == 0.0)
    {
      fprintf(stderr,"Value of Alpha resolution ( = 0.0) not allowed.\n");
      exit(-1);
    }      

  /* don't allow negative bands (for safty in griding-routines) */
  if ( (uvar_AlphaBand < 0) ||  (uvar_DeltaBand < 0) )
    {
      fprintf (stderr, "\nNegative value of sky-bands not allowed (alpha or delta)!\n");
      exit(-1);
    }

  /* check consistency of bounds for the sky-region */
  /*  FIXME: actually that's nonsense, as the region might stretch across this coordinate "cut" */
  if ( (uvar_Alpha < 0) || (uvar_Alpha >= LAL_TWOPI) || (fabs(uvar_Delta) > LAL_PI_2) 
       ||(uvar_Alpha + uvar_AlphaBand >= LAL_TWOPI) || (uvar_Delta + uvar_DeltaBand > LAL_PI_2) )
    {
      fprintf (stderr, "\nSorry, sky-region is not supposed to cross either 0-meridian or a pole!\n");
      exit(-1);
    }


  strcpy(cfg->EphemEarth,uvar_EphemDir);
  strcat(cfg->EphemEarth,"/earth");
  sprintf (tmp, "%02d", uvar_EphemYear);
  strcat(cfg->EphemEarth, tmp);
  strcat(cfg->EphemEarth,".dat");

  strcpy(cfg->EphemSun,uvar_EphemDir);
  strcat(cfg->EphemSun,"/sun");
  strcat(cfg->EphemSun, tmp);
  strcat(cfg->EphemSun,".dat");

  /* *** Make sure the e-files are really there *** */
      fp=fopen(cfg->EphemEarth,"r");
      if (fp==NULL) 
	{
	  fprintf(stderr,"Could not find %s\n",cfg->EphemEarth);
	  return 1;
	}
      fclose(fp);
      fp=fopen(cfg->EphemSun,"r");
      if (fp==NULL) 
	{
	  fprintf(stderr,"Could not find %s\n",cfg->EphemSun);
	  return 1;
	}
      fclose(fp);
  /* ********************************************** */

  strcpy(command,uvar_DataDir);
  strcat(command,"/*");
  /* this default-value needs to be set-up here */
  if (uvar_BaseName == NULL) {
    uvar_BaseName = LALMalloc (10);
    strcpy (uvar_BaseName, "SFT");
  }
  strcat(command,uvar_BaseName);
  strcat(command,"*");
  
  globbuf.gl_offs = 1;
  glob(command, GLOB_ERR, NULL, &globbuf);

  /* read file names -- MUST NOT FORGET TO PUT ERROR CHECKING IN HERE !!!! */

  if(globbuf.gl_pathc==0)
    {
      fprintf(stderr,"No SFTs in directory %s ... Exiting.\n",uvar_DataDir);
      return 1;
    }

  while ((UINT4)fileno < globbuf.gl_pathc) 
    {
      strcpy(cfg->filelist[fileno],globbuf.gl_pathv[fileno]);
      fileno++;
      if (fileno > MAXFILES)
	{
	  fprintf(stderr,"Too many files in directory! Exiting... \n");
	  return 1;
	}
    }
  globfree(&globbuf);

  cfg->SFTno=fileno; /* remember this is 1 more than the index value */

  /* initialize detector */
  if(uvar_IFO == 0) cfg->Detector=lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  if(uvar_IFO == 1) cfg->Detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
  if(uvar_IFO == 2) cfg->Detector=lalCachedDetectors[LALDetectorIndexLHODIFF];
  if(uvar_IFO == 3)
    {
        if (CreateDetector(&(cfg->Detector))) return 5;
    }


  /* open FIRST file and get info from it*/
  fp=fopen(cfg->filelist[0],"r");
  /* read in the header from the file */
  errorcode=fread((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1) 
    {
      fprintf(stderr,"No header in data file %s\n",cfg->filelist[0]);
      return 1;
    }

  /* check that data is correct endian order */
  if (header.endian!=1.0)
    {
      fprintf(stderr,"First object in file %s is not (double)1.0!\n",cfg->filelist[0]);
      fprintf(stderr,"It could be a file format error (big/little\n");
      fprintf(stderr,"endian) or the file might be corrupted\n\n");
      return 2;
    }
  fclose(fp);

  cfg->Ti=header.gps_sec;  /* INITIAL TIME */

  /* open LAST file and get info from it*/
  fp=fopen(cfg->filelist[fileno-1],"r");
  /* read in the header from the file */
  errorcode=fread((void*)&header,sizeof(header),1,fp);
  if (errorcode!=1) 
    {
      fprintf(stderr,"No header in data file %s\n",cfg->filelist[fileno-1]);
      return 1;
    }
  /* check that data is correct endian order */
  if (header.endian!=1.0)
    {
      fprintf(stderr,"First object in file %s is not (double)1.0!\n",cfg->filelist[fileno-1]);
      fprintf(stderr,"It could be a file format error (big/little\n");
      fprintf(stderr,"endian) or the file might be corrupted\n\n");
      return 2;
    }
  fclose(fp);

  cfg->Tf=header.gps_sec+header.tbase;  /* FINAL TIME */

  cfg->tsft=header.tbase;  /* Time baseline of SFTs */
    
  /* if user has not input demodulation frequency resolution; set to 1/Tobs */
  if( uvar_dFreq == 0.0 ) 
    uvar_dFreq=1.0/(2.0*header.tbase*cfg->SFTno);

  cfg->FreqImax=(INT4)(uvar_FreqBand/uvar_dFreq+.5)+1;  /*Number of frequency values to calculate F for */
    
  /* if user has not input demodulation frequency resolution; set to 1/Tobs */
  if( uvar_dSpin == 0.0 ) uvar_dSpin=1.0/(2.0*header.tbase*cfg->SFTno*(cfg->Tf-cfg->Ti));

  cfg->SpinImax=(int)(uvar_SpinBand/uvar_dSpin+.5)+1;  /*Number of frequency values to calculate F for */

  cfg->nsamples=header.nsamples;    /* # of freq. bins */

  /* frequency resolution: used only for printing! */
  df=(1.0)/(1.0*header.tbase);
  cfg->df=df;

  cfg->ifmax=ceil((1.0+DOPPLERMAX)*(uvar_Freq+uvar_FreqBand)*cfg->tsft)+uvar_Dterms;
  cfg->ifmin=floor((1.0-DOPPLERMAX)*uvar_Freq*cfg->tsft)-uvar_Dterms;

  /* allocate F-statistic arrays */
  Fstat.F =(REAL8*)LALMalloc(cfg->FreqImax*sizeof(REAL8));
  if(uvar_EstimSigParam) 
    {
      Fstat.Fa =(COMPLEX16*)LALMalloc(cfg->FreqImax*sizeof(COMPLEX16));
      Fstat.Fb =(COMPLEX16*)LALMalloc(cfg->FreqImax*sizeof(COMPLEX16));
    } else {
      Fstat.Fa = NULL;
      Fstat.Fb = NULL;
    }

  /* Tell the user what we have arrived at */
#ifdef DEBG_SGV
    fprintf(stdout,"\n");
    fprintf(stdout,"# SFT time baseline:                  %f min\n",header.tbase/60.0);
    fprintf(stdout,"# SFT freq resolution:                %f Hz\n",df);
    fprintf(stdout,"# Starting search frequency:          %f Hz\n",uvar_Freq);
    fprintf(stdout,"# Demodulation frequency band:        %f Hz\n",uvar_FreqBand);
    fprintf(stdout,"# no of SFT in a DeFT:                %f\n",ceil((1.0*(cfg->Tf - cfg->Ti))/header.tbase));
    fprintf(stdout,"# Actual # of SFTs:                   %d\n",cfg->SFTno);
    fprintf(stdout,"# ==> DeFT baseline:                  %f hours\n",(cfg->Tf - cfg->Ti)/3600.0);
#endif

  return 0;  
}

/*******************************************************************************/



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

  /* Free SFTData */
  for (k=0;k<GV.SFTno;k++)
    {
      LALFree(SFTData[k]->fft->data->data);
      LALFree(SFTData[k]->fft->data);
      LALFree(SFTData[k]->fft);
      LALFree(SFTData[k]);
    }
  LALFree(SFTData);

  /* Free timestamps */
  LALFree(timestamps);

  LALFree(Fstat.F);
  if(uvar_EstimSigParam) 
    {
      LALFree(Fstat.Fa);
      LALFree(Fstat.Fb);
    }

  /* Free DemodParams */
  LALSDestroyVector(&status, &(DemodParams->amcoe->a));
  LALSDestroyVector(&status, &(DemodParams->amcoe->b));

  LALFree(DemodParams->skyConst);
  LALFree(DemodParams->spinDwn);
  LALFree(DemodParams);

  /* Free config-Variables and userInput stuff */
  LALFreeUserVars (&status, uvars);

  /* Free DopplerScan-stuff (grid) */
  FreeDopplerScan (&status, &thisScan);

  /* this comes from clusters.c */
  if (highFLines->clusters) LALFree(highFLines->clusters);
  if (highFLines->Iclust) LALFree(highFLines->Iclust);
  if (highFLines->NclustPoints) LALFree(highFLines->NclustPoints);

  /* did we forget anything ? */
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
		  uvar_Freq+indexes[ntop]*uvar_dFreq,
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
  err=fprintf(fpmax,"%10.5f %10.8f %10.8f    %d %10.5f %10.5f %10.5f\n",uvar_Freq,
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
  
  REAL4 xre,xim;

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
       REAL4 freq;
       REAL8 r0,r1;
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
     REAL4 freq;
     REAL8 r0,r1,r2;
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
     REAL4 freq;
     REAL8 r0,r1;
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
  /* INT2 windowSize=(0.01/uvar_dFreq);               0.1 is 1E-4*1000 */
  INT2 windowSize=100;
  REAL4 THR=10.0;
  
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

  THR=uvar_Fthreshold;


  /* wings=windowSize/2; */
  /*  0.0002 is the max expected width of the F stat curve for signal */
  /*  with ~ 10 h observation time */
  /*  0.0001 = 0.0002/2 */
  /*  let me put 0.005 */
  dmp=0.5+0.0002/uvar_dFreq;
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
  outliersInput->ifmin=((uvar_Freq/uvar_dFreq)+0.5);
  outliersInput->data = F1;

  ComputeOutliers(outliersInput, outliersParams, outliers);

   if (outliers->Noutliers == 0){

#ifdef FILE_FTXT
     /*  F.txt file contains freq, F, noise floor of F   */
     for (i=0;i<nbins;i++){ 
       REAL4 freq;
       REAL8 r0,r1;
       freq=uvar_Freq + i*uvar_dFreq;
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
     REAL4 freq;
     REAL8 r0,r1,r2;
     j=SpLines->Iclust[i];
     freq=(uvar_Freq+SpLines->Iclust[i]*uvar_dFreq);
     r0=F1->data[j];
     r1=FloorF1->data[j];
     r2=SpLines->clusters[i]*FloorF1->data[j];
     fprintf(outfile1,"%f %E %E %E\n",freq,r0,r1,r2);
   }
#endif
#ifdef FILE_FTXT   
   /*  PSD.txt file contains freq, PSD, noise floor   */
   for (i=0;i<nbins;i++){ 
     REAL4 freq;
     REAL8 r0,r1;
     freq=uvar_Freq + i*uvar_dFreq;
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

  if (uvar_SignalOnly != 1)
    RngMedBias (&status, &medianbias, windowSize);


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
      
      if(uvar_SignalOnly == 1){
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

