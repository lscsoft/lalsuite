/*<lalVerbatim file="LALInspiralParseParametersCV"> 
  Author:  Cokelaer, T.
  $Id$
  </lalVerbatim>  */



/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralParseParameters.c}}
DOC IN PROGRESS


Module to work with the inspiralTemplate Structure. 

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralITStructureParseParametersCP}
\input{LALInspiralITStructurePrintCP}
\input{LALInspiralITStructureSetDefaultCP}
\input{LALInspiralITStructureInit2DummyCP}
\input{LALInspiralITStructureHelpCP}
\input{LALInspiralITStructureCheckCP}


\idx{LALInspiralITStructureParseParameters}
\idx{LALInspiralITStructurePrint}
\idx{LALInspiralITStructureSetDefault}
\idx{LALInspiralITStructureInit2Dummy}
\idx{LALInspiralITStructureHelp}
\idx{LALInspiralITStructureCheck}


\idx{LALInspiralParseParametersForInjection}




\subsubsection*{Description}
This module is a set of functions to play with the inspiralTemplate structure of the 
inspiral package.

\begin{itemize}
\item  The \texttt{LALInspiralITSructureParseParameters} function allows 
the user to parsea string with respect to that structure. Each variable in the 
inspiralTemplate structure might be parse with a string like "--(name of the variable)+(value)"
 i.e. \textit{--approximant TaylorT1}. Each argument starts with a double dash character
followed by a key word which is exactly as written in the InspiralTemplate Structure such as
 --order, --mass1, --mass2, --fCutoff ...

Once the string is parsed, the checking function is called. 
 
\item The \texttt{LALInspiralITStructurePrint} function will print on the stdout the value 
of the InspiralTemplate structure.
\item The \texttt{LALInspiralITStructureSetDefault} set default values to the variables. 
Those values are written in the C-code.
\item \texttt{LALInspiralITStructureInit2Dummy}
\item \texttt{LALInspiralITStructureCheck}
\item \texttt{LALInspiralITStructureHelp}
\end{itemize}


\subsubsection*{Algorithm}
None

\vfill{\footnotesize\input{LALInspiralParseParametersCV}}

</lalLaTeX>  */

#include <lal/LALInspiral.h>
#include <lal/StringInput.h>



#define  INSPIRALTEMPLATE_APPROXIMANT 	TaylorT3
#define  INSPIRALTEMPLATE_ORDER 	twoPN
#define  INSPIRALTEMPLATE_MASS1  	10.
#define  INSPIRALTEMPLATE_MASS2 	10.
#define  INSPIRALTEMPLATE_FCUTOFF 	1000.
#define  INSPIRALTEMPLATE_FLOWER 	40.
#define  INSPIRALTEMPLATE_TSAMPLING 	2048.
#define  INSPIRALTEMPLATE_DISTANCE 	1.   /*MPC*/
#define  INSPIRALTEMPLATE_SIGNALAMPLITUDE 1.
#define  INSPIRALTEMPLATE_STARTPHASE    0.
#define  INSPIRALTEMPLATE_STARTTIME     0.

#define  INSPIRALTEMPLATE_THETA 	0.
#define  INSPIRALTEMPLATE_ZETA2 	0.
#define  INSPIRALTEMPLATE_OMEGAS	0.

#define  INSPIRALTEMPLATE_ALPHA 	0.
#define  INSPIRALTEMPLATE_PSI0 		100000.
#define  INSPIRALTEMPLATE_PSI3 		-1000.

#define  INSPIRALTEMPLATE_ALPHA1 	0.
#define  INSPIRALTEMPLATE_ALPHA2 	0.
#define  INSPIRALTEMPLATE_BETA 		0.

#define  INSPIRALTEMPLATE_INCLINATION 	0.
#define  INSPIRALTEMPLATE_ECCENTRICITY 	0.
#define  INSPIRALTEMPLATE_ORBITTHETA0 	0.
#define  INSPIRALTEMPLATE_ORBITPHI0 	0.
#define  INSPIRALTEMPLATE_SOURCETHETA 	0.
#define  INSPIRALTEMPLATE_SOURCEPHI 	0.

NRCSID (LALINSPIRALPARSEPARAMETERSC, "");


void LALInspiralITStrucutreInit2Dummy(LALStatus *status, 
				      InspiralTemplate *params);



/* <lalVerbatim file="LALInspiralITStructureParseParametersCP"> */
void LALInspiralITStructureParseParameters(LALStatus *status,
					   UINT4 argc, 
					   CHAR **argv,
					   InspiralTemplate *params)
/* </lalVerbatim> */
{
  UINT4 i	= 1;


  INITSTATUS( status, "LALInspiralParseParameters", LALINSPIRALPARSEPARAMETERSC);
  ATTATCHSTATUSPTR( status );
  
  while(i <argc)
    {
      if (strcmp(argv[i], "--approximant")==0)
	{
	  params->massChoice  = m1Andm2;
	  if (strcmp(argv[++i],"TaylorT1")==0){
	    params->approximant = TaylorT1; }
	  else if (strcmp(argv[i],"TaylorT2")==0){
	    params->approximant = TaylorT2;}
	  else if (strcmp(argv[i],"TaylorT3")==0){
	    params->approximant = TaylorT3;}
	  else if (strcmp(argv[i],"TaylorF1")==0){
	    params->approximant = TaylorF1;}
	  else if (strcmp(argv[i],"TaylorF2")==0){
	    params->approximant = TaylorF2;}
	  else if (strcmp(argv[i],"PadeT1")==0){
	    params->approximant = PadeT1;}
	  else if (strcmp(argv[i],"PadeF1")==0){
	    params->approximant = PadeF1;}
	  else if (strcmp(argv[i],"EOB")==0){
	    params->approximant = EOB;}
	  else if (strcmp(argv[i],"BCV")==0){
	    params->approximant = BCV;
	    params->massChoice  = psi0Andpsi3;
	  }
	  else if (strcmp(argv[i],"BCVSpin")==0){
	    params->massChoice  = psi0Andpsi3;
	    params->approximant = BCVSpin;
	  }
	  else if (strcmp(argv[i],"SpinTaylorT3")==0){
	    params->approximant = SpinTaylorT3; }
	  else {params->approximant = TaylorT1;} /*is it correct ? */
	}/* SpinTaylor is not available here only for inject package*/
      else if (strcmp(argv[i], "--order")==0){
	params->order = atoi(argv[++i]);}
      else if (strcmp(argv[i],"--mass1")==0){
	params->mass1 = atof(argv[++i]);}
      else if (strcmp(argv[i],"--mass2")==0){
	params->mass2 = atof(argv[++i]);}
      else if (strcmp(argv[i],"--massChoice")==0){
	if (strcmp(argv[++i],"masses")==0){
	  params->massChoice = m1Andm2;}
	else if (strcmp(argv[i],"psi")==0){ 
	  params->massChoice = psi0Andpsi3;}
      }
      else if (strcmp(argv[i],"--fCutoff")==0){
	params->fCutoff = atof(argv[++i]);}
      else if (strcmp(argv[i],"--fFinal")==0){
	params->fFinal = atof(argv[++i]);}
      else if (strcmp(argv[i],"--fLower")==0){
	params->fLower = atof(argv[++i]);}
      else if (strcmp(argv[i],"--tSampling")==0){
	params->tSampling = atof(argv[++i]); }
      else if (strcmp(argv[i],"--distance")==0){
	params->distance = atof(argv[++i]); }
      else if (strcmp(argv[i],"--startPhase")==0){
	params->signalAmplitude = atof(argv[++i]); }
      else if (strcmp(argv[i],"--startTime")==0){
	params->startTime = atof(argv[++i]); }
      else if (strcmp(argv[i],"--Theta")==0){
	params->Theta = atof(argv[++i]); }
      else if (strcmp(argv[i],"--Zeta2")==0){
	params->Zeta2 = atof(argv[++i]); } 
      else if (strcmp(argv[i],"--alpha")==0){
	params->alpha = atof(argv[++i]); }
      else if (strcmp(argv[i],"--alpha1")==0){
	params->alpha1 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--alpha2")==0){
	params->alpha2 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--beta")==0){
	params->beta = atof(argv[++i]); }  
      else if (strcmp(argv[i],"--psi0")==0){
	params->psi0 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--psi3")==0){
	params->psi3 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--inclination")==0){
	params->inclination = atof(argv[++i]); }
      else if (strcmp(argv[i],"--orbitTheta0")==0){
	params->orbitTheta0 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--orbitPhi0")==0){
	params->orbitPhi0 = atof(argv[++i]); }
      else if (strcmp(argv[i],"--sourceTheta")==0){
	params->sourceTheta = atof(argv[++i]); }
      else if (strcmp(argv[i],"--sourcePhi")==0){
	params->sourcePhi = atof(argv[++i]); }
      /* we need also the spin .Has to be checked*/
      else if (strcmp(argv[i],"--spin")==0)
	{	 
	  params->spin1[0]= atof(argv[++i]);
	  params->spin2[0]= atof(argv[++i]);
	}
      else if (strcmp(argv[i],"--angle")==0)
	{	 
	  params->spin1[1]= atof(argv[++i]);
	  params->spin1[2]= atof(argv[++i]);
	  params->spin2[1]= atof(argv[++i]);
	  params->spin2[2]= atof(argv[++i]);
	  
	}
      else if (strcmp(argv[i],"--eccentricity")==0){
	params->sourcePhi = atof(argv[++i]); }
       /* if the option is not recongnized, we skip the next 
		   argument and go ahead. We don't want to quit */      
      i++;       
   }
  
  LALInspiralITStructureCheck(status->statusPtr, *params);
  CHECKSTATUSPTR(status);
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
}


/* <lalVerbatim file="LALInspiralITStructureInit2DummyCP"> */
void LALInspiralITStrucutreInit2Dummy(LALStatus *status, 
				      InspiralTemplate *params)
/* </lalVerbatim> */
{
  params->approximant  = -1;
  params->order        = -1;
  params->mass1        = -1;
  params->mass2        = -1;
  params->fCutoff      = -1;
  params->fLower       = -1;
  params->tSampling    = -1; 
  params->distance     = -1; 
  params->signalAmplitude = -1;
  params->startPhase   =-1;
  params->startTime    =-1;
  
  params->Theta        = -1;
  params->Zeta2        = -1;

  params->alpha        = -1;
  params->psi0         = -1;
  params->psi3         = -1;
  
  params->alpha1       = -1;
  params->alpha2       = -1;
  params->beta         = -1;
  
  params->Theta        = -1;
  params->Zeta2        = -1;

  params->inclination  = -1;
  params->orbitTheta0  = -1;
  params->orbitPhi0    = -1;
  params->sourceTheta  = -1;
  params->sourcePhi    = -1;
  
  params->eccentricity  = -1;  
  RETURN(status);

}


/* <lalVerbatim file="LALInspiralITStructureCheckCP"> */
void LALInspiralITStructureCheck(LALStatus *status, 
				 InspiralTemplate params)
/* </lalVerbatim> */
{

  INITSTATUS( status, "LALInspiralITSStrucutreCheck", LALINSPIRALPARSEPARAMETERSC);
  ATTATCHSTATUSPTR( status );

  switch (params.approximant)
    {
    case BCV:
/*      if (params.alpha < 0) {  
	ABORT(status, LALINSPIRALH_EALPHA, LALINSPIRALH_MSGEALPHA);
      }*/
       if (params.fCutoff < 0){
	ABORT(status, LALINSPIRALH_EFCUTOFF,  LALINSPIRALH_MSGEFCUTOFF);
      }
      else if (params.psi0 < 0) {
	ABORT(status, LALINSPIRALH_EPSI0,  LALINSPIRALH_MSGEPSI0);
      }
      else if (params.psi3 > 0) {
	ABORT(status, LALINSPIRALH_EPSI3,  LALINSPIRALH_MSGEPSI3);
      }
       
      break;
    case PadeF1:    
      ABORT(status, LALINSPIRALH_EAPPROXIMANT, LALINSPIRALH_MSGEAPPROXIMANT );
      break;
    case TaylorT1: 
    case TaylorT2: 
    case TaylorT3:
    case BCVSpin:
    case PadeT1:
    case TaylorF1:
    case EOB:
    case TaylorF2:
    case SpinTaylorT3:
    case SpinTaylor:
      break;
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);    
}

/* <lalVerbatim file="LALInspiralITStructurePrintCP"> */
void LALInspiralITStructurePrint(LALStatus *status, 
				 InspiralTemplate params)
/* </lalVerbatim> */
{
  /* later */
  printf("# approximant = %-15.12d\n", params.approximant);
  printf("# order       = %-15.12d\n", params.order);
  printf("# mass1       = %-15.12f\n", params.mass1); 
  printf("# mass2       = %-15.12f\n", params.mass2);
  printf("# fcutoff     = %-15.12f\n", params.fCutoff);
  printf("# fLower      = %-15.12f\n", params.fLower);
  printf("# tSampling   = %-15.12f\n", params.tSampling);
  printf("# distance    = %-15.12f\n", params.distance);
  printf("# startPhase  = %-15.12f\n", params.startPhase);
  printf("# startTime   = %-15.12f\n", params.startTime);

  printf("# zeta2       = %-15.12f\n", params.Zeta2);
  printf("# omegaS      = %-15.12f\n", params.OmegaS);

  printf("# alpha       = %-15.12f\n", params.alpha);
  printf("# psi0        = %-15.12f\n", params.psi0);
  printf("# psi3        = %-15.12f\n", params.psi3);
  printf("# alpha1      = %-15.12f\n", params.alpha1);
  printf("# alpha2      = %-15.12f\n", params.alpha2);
  printf("# beta        = %-15.12f\n", params.beta);

  printf("# inclination = %-15.12f\n", params.inclination);
  printf("# orbitTheta0 = %-15.12f\n", params.orbitTheta0);
  printf("# orbitPhi0   = %-15.12f\n", params.orbitPhi0);
  printf("# sourceTheta = %-15.12f\n", params.sourceTheta);
  printf("# sourcePhi   = %-15.12f\n", params.sourcePhi);

  printf("# eccentricity= %-15.12f\n", params.eccentricity);



/* Paramters which are computed using LALInspiralParameterCalc */

  printf("# chirpMass   = %-15.12f\n", params.chirpMass); 
  printf("# eta         = %-15.12f\n", params.eta);
  printf("# totalMass   = %-15.12f\n", params.totalMass); 
  printf("# fFinal      = %-15.12f\n", params.fFinal);
  printf("# t0          = %-15.12f\n", params.t0); 
  printf("# t2          = %-15.12f\n", params.t2); 
  printf("# t3          = %-15.12f\n", params.t3); 
  printf("# t4          = %-15.12f\n", params.t4); 
  printf("# t5          = %-15.12f\n", params.t5); 
  printf("# tC          = %-15.12f\n", params.tC); 

  printf("# massChoice  = %-15.12d\n", params.massChoice); 

  RETURN(status);
}


/**/
/* <lalVerbatim file="LALInspiralITStructureSetDefaultCP"> */
void LALInspiralITStructureSetDefault(LALStatus *status,
				      InspiralTemplate *params)
/* </lalVerbatim> */
{
  
  params->approximant  = INSPIRALTEMPLATE_APPROXIMANT;
  params->order        = INSPIRALTEMPLATE_ORDER ;
  params->mass1        = INSPIRALTEMPLATE_MASS1 ;
  params->mass2        = INSPIRALTEMPLATE_MASS2 ;
  params->fCutoff      = INSPIRALTEMPLATE_FCUTOFF;
  params->fLower       = INSPIRALTEMPLATE_FLOWER;
  params->tSampling    = INSPIRALTEMPLATE_TSAMPLING; 
  params->distance     = INSPIRALTEMPLATE_DISTANCE;   
  
  params->signalAmplitude = INSPIRALTEMPLATE_SIGNALAMPLITUDE;
  params->startPhase   = INSPIRALTEMPLATE_STARTPHASE;
  params->startTime    = INSPIRALTEMPLATE_STARTTIME;
  
  params->Theta        = INSPIRALTEMPLATE_THETA;
  params->Zeta2        = INSPIRALTEMPLATE_ZETA2;  
  
  
  params->alpha        = INSPIRALTEMPLATE_ALPHA;
  params->psi0         = INSPIRALTEMPLATE_PSI0;
  params->psi3         = INSPIRALTEMPLATE_PSI3;
   

  params->alpha1       = INSPIRALTEMPLATE_ALPHA1;
  params->alpha2       = INSPIRALTEMPLATE_ALPHA2;
  params->beta         = INSPIRALTEMPLATE_BETA;
    

  params->inclination  = INSPIRALTEMPLATE_INCLINATION;
  params->orbitTheta0  = INSPIRALTEMPLATE_ORBITTHETA0;
  params->orbitPhi0    = INSPIRALTEMPLATE_ORBITPHI0;
  params->sourceTheta  = INSPIRALTEMPLATE_SOURCETHETA;
  params->sourcePhi    = INSPIRALTEMPLATE_SOURCEPHI;
    
  params->eccentricity = INSPIRALTEMPLATE_ECCENTRICITY;   
  
  params->ieta         = 1.;  
  params->OmegaS       = INSPIRALTEMPLATE_OMEGAS;
  params->nStartPad    = 0;
  params->nEndPad      = 0;
  params->massChoice   = m1Andm2;  

  RETURN(status);
}


/* <lalVerbatim file="LALInspiralITStructureHelpCP"> */
void LALInspiralITStructureHelp()
/* </lalVerbatim> */
{

  fprintf(stderr,"InspiralTemplate Structure; parsing arguments\n");
  fprintf(stderr,"--approximant (TaylorT1, TaylorT2, TaylorT3, EOB, BCV, BCVSpin, PadeT1)\n");
  fprintf(stderr,"--order       (0, 1, 2, 3, 4, 5, 6 (i.e. 4==twoPN)\n");
  fprintf(stderr,"--mass1       (in solar mass)\n"); 
  fprintf(stderr,"--mass2       (in solar mass)\n");
  fprintf(stderr,"--fcutoff     \n");
  fprintf(stderr,"--fLower      \n");
  fprintf(stderr,"--tSampling   \n");
  fprintf(stderr,"--distance    \n");
  fprintf(stderr,"--startPhase  \n");
  fprintf(stderr,"--startTime   \n");
  fprintf(stderr,"--zeta2       (EOB 3PN)\n");
  fprintf(stderr,"--omegaS      (EOB 3PN)\n");

  fprintf(stderr,"--alpha       (BCV must be > 0)\n");
  fprintf(stderr,"--psi0        (BCV must be > 0)\n");
  fprintf(stderr,"--psi3        (BCV must be < 0)\n");
  fprintf(stderr,"--alpha1      (BCVSPin)\n");
  fprintf(stderr,"--alpha2      (BCVSPin)\n");
  fprintf(stderr,"--beta        (BCVSPIN)\n");

  fprintf(stderr,"--inclination \n");
  fprintf(stderr,"--orbitTheta0 \n");
  fprintf(stderr,"--orbitPhi0   \n");
  fprintf(stderr,"--sourceTheta \n");
  fprintf(stderr,"--sourcePhi   \n");

  fprintf(stderr,"--eccentricity\n");

}
