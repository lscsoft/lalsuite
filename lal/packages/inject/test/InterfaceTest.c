/***************************************** <lalVerbatim file="InterfaceTest">
Author: Cokelaer, T.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>


DRAFT DOCUMENTATION
 right now this test file read an input file and for each valid line compute
 the approriate waveform. The waveform are eiter produce inside the inject pacakge 
either inside the inspiral package. 

This file is composed of two main part. One to read the file and parse the parameter and one
to generate the waveform.

Parsing the parameters has required to create a new struct GeneralInspiralStruc which can parse either
"inspiral" data either "inject" data depending of the input data. 

Then the function LALGeneralInspiral create the amplitude, freq and phi vectors needed for further
injection including the data itself (noise, h(t)) 



\subsection{Program \texttt{InterfaceTest.c}}
\label{ss:InterfaceTest.c}

Interface to generate any kind of gravitational waves signal.

\subsubsection*{Usage}
\begin{verbatim}
InterfaceTest [-p parameterfile] 
\end{verbatim}

\subsubsection*{Description}


\subsubsection*{Exit codes}
****************************************** </lalLaTeX><lalErrTable> */



#define INTERFACETESTC_ENORM  0
#define INTERFACETESTC_ESUB   1
#define INTERFACETESTC_EARG   2
#define INTERFACETESTC_EVAL   3
#define INTERFACETESTC_EFILE  4
#define INTERFACETESTC_EINPUT 5
#define INTERFACETESTC_EMEM   6

#define INTERFACETESTC_MSGENORM  "Normal exit"
#define INTERFACETESTC_MSGESUB   "Subroutine failed"
#define INTERFACETESTC_MSGEARG   "Error parsing arguments"
#define INTERFACETESTC_MSGEVAL   "Input arg"
#define INTERFACETESTC_MSGEFILE  "Could not open file"
#define INTERFACETESTC_MSGEINPUT "Error reading file"
#define INTERFACETESTC_MSGEMEM   "Out of memory"

/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
liste des fonctions utilisees
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{InterfaceTestCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALInspiral.h>
#include <lal/AVFactories.h>
#include <lal/StringInput.h>
#include <lal/SeqFactories.h>
#include <lal/StreamInput.h>
#include <lal/StreamOutput.h>
#include <lal/LALInitBarycenter.h>
#include <lal/FrameData.h>
#include <lal/LALConfig.h>
#include <lal/FrameStream.h>


NRCSID( INTERFACETESTC, "$Id$" );

/* Default value for the generation of GW waveforms*/
#define INTERFACETEST_APPROXIMANT      PPN
#define INTERFACETEST_M1          1.4
#define INTERFACETEST_M2         1.4
#define INTERFACETEST_DISTANCE   10  /*mpc*/
#define INTERFACETEST_TC          0.
#define INTERFACETEST_ORDER       4
#define INTERFACETEST_INCLINATION 0.
#define INTERFACETEST_DEC         0.
#define INTERFACETEST_RA          0.
#define INTERFACETEST_PSI         0.
#define INTERFACETEST_PHIC        0.
#define INTERFACETEST_FLOWER     30.
#define INTERFACETEST_FUPPER   ((INTERFACETEST_SAMPLING-1)/2.)   /*must be < sampling/2*/
#define INTERFACETEST_ZETA2       0.
#define INTERFACETEST_OMEGAS      0.
#define INTERFACETEST_THETA       0.
#define INTERFACETEST_T0          0.
#define INTERFACETEST_T1          1.
#define INTERFACETEST_T2          2.
#define INTERFACETEST_SAMPLING 2048.
#define INTERFACETEST_A1          1000.
#define INTERFACETEST_A2          1000.
#define INTERFACETEST_PHI0        0.
#define INTERFACETEST_F0         100.
#define INTERFACETEST_ARG         0
#define INTERFACETEST_UDOT        0.5
#define INTERFACETEST_RP          0.5
#define INTERFACETEST_E           0.   
#define INTERFACETEST_ALPHA       0.
#define INTERFACETEST_FBCV        INTERFACETEST_FUPPER
#define INTERFACETEST_DELTAT      0  /*default use deltaT = 1/sampling*/
#define INTERFACETEST_LENGTH   1024  
#define INTERFACETEST_LENGTHIN  1e9

                                                                                                                                                                            


#define WHITESPACE " \f\n\r\t\v" /* list of whitespace characters */

#define MSGLEN 1024  /* Max. length of warning/info messages */

/* Usage format string. */
#define USAGE "Usage: %s [-p paramfile] \n"


static void ClearStatus (LALStatus *status);
void ReadLine(LALStatus *status, 
	      FILE *fp,
	      int *argc,
	      GeneralInspiralStruc *);
int  ParametersParse(LALStatus *, 
		     int , 
		     char **, 
		     GeneralInspiralStruc *);
void Help();

INT4 lalDebugLevel=1;



int main(int argc, char **argv)
{
  GeneralInspiralStruc params;
  CoherentGW   waveform ; 
  static   LALStatus status;
  
  
  CHAR *paramfile  = NULL;
  INT4 i=1;
  INT4 lineargc=0  ;
  
  FILE *fp =NULL;;
 
  

  
  if (argc!=3)
    {
      return INTERFACETESTC_EARG;
    }
  else
    {
      while( i<argc)
	{
	  if (strcmp(argv[i], "-p")==0)
	    paramfile = argv[++i];
	  i++;
	}
      
    }
  
  
  if (paramfile){ 
    if ( ( fp = fopen( paramfile, "r" ) ) == NULL ) {   
      return INTERFACETESTC_EFILE;
    }
  }
  
  
/*    while the file contains a new line we extract the parameters and try to inject the signal */
/*      into ~some~ data.    */
  
  while( !feof(fp) ) 
    {        
      memset(&params,0,sizeof(GeneralInspiralStruc));
      ReadLine(&status, fp, &lineargc, &params); 
      ClearStatus(&status);
      printf("#number of parameters %d\n", lineargc);      

       if (lineargc >0) 
 	{	 
	  memset(&waveform,0,sizeof(CoherentGW));
	  LALGenerateInspiral(&status, &waveform, &params);   
	  ClearStatus(&status);
	  /*print */
	  for (i=0; i<waveform.f->data->length; i++) 
	  printf("%e\n",waveform.f->data->data[i]); 
 	  printf("&\n"); 


	  /*free memory*/
	  if (params.inspiral.approximant !=TaylorT2 || params.inspiral.approximant!=TaylorT3)
	    {
	      LALSDestroyVectorSequence( &status, &(waveform.a->data)); 
	      ClearStatus(&status);
	      LALSDestroyVector( &status, &(waveform.f->data) );
	      ClearStatus(&status);
	      LALDDestroyVector( &status, &(waveform.phi->data)); 
	      ClearStatus(&status);
	      
	      LALFree( waveform.a );
	      LALFree( waveform.f );
	      LALFree( waveform.phi );
	    }
	  
	  if (params.inspiral.approximant==SpinOrbitCW)
	  {
	    LALDDestroyVector( &status, &(params.socw.f));
	    ClearStatus(&status);
	  }
	  if (params.inspiral.approximant==TaylorCW)
	    {
	      LALDDestroyVector( &status, &(params.taylorcw.f));
	      ClearStatus(&status);
	    }
	  
	  
	  
	}
       lineargc = 0;              
    } 

  
  
  
  /*todo check leak memory and co*/
  
  
  LALCheckMemoryLeaks();

  return 0;
}





void ReadLine(
	      LALStatus *status,
	      FILE *fp, 
	      int *argc,
	      GeneralInspiralStruc *params)
{

  CHARVector *line = NULL ;
  TokenList *tokens = NULL; /* input line parsed into tokens */


  LALCHARReadVector( status, &line, fp );
  ClearStatus (status);

    long i;
  for ( i = 0; i < line->length; i++ )
    if ( line->data[i] == '%' || line->data[i] == '#' ) {
      line->data[i] = '\0';
      i = line->length;
    }
  
  LALCreateTokenList( status, &tokens, line->data, WHITESPACE);
  ClearStatus (status);
  LALCHARDestroyVector( status, &line);
  ClearStatus (status);
  
  
  *argc = tokens->nTokens;
  if (tokens->nTokens >0)
    {      
      ParametersParse(status,
		      tokens->nTokens,
		      tokens->tokens,
		      params);
    }
  
  LALDestroyTokenList(status, &tokens);
  ClearStatus (status);
   
  
}

int  ParametersParse(LALStatus *status,
		     int   argc, 
		     char **argv,
		     GeneralInspiralStruc *params)
{
  REAL8 
    tc       = INTERFACETEST_TC,
    m1       = INTERFACETEST_M1,
    m2       = INTERFACETEST_M2,
    d        = INTERFACETEST_DISTANCE,
    inc      = INTERFACETEST_INCLINATION,
    ra       = INTERFACETEST_RA,
    dec      = INTERFACETEST_DEC,
    psi      = INTERFACETEST_PSI,
    phic     = INTERFACETEST_PHIC,
    fi       = INTERFACETEST_FLOWER,
    ff       = INTERFACETEST_FUPPER,
    zeta2    = INTERFACETEST_ZETA2,
    omegas   = INTERFACETEST_OMEGAS,
    theta    = INTERFACETEST_THETA,
    t0       = INTERFACETEST_T0,
    t1       = INTERFACETEST_T1,
    t2       = INTERFACETEST_T2,
    a1       = INTERFACETEST_A1,
    a2       = INTERFACETEST_A2,
    phi0     = INTERFACETEST_PHI0,
    f0       = INTERFACETEST_F0,
    arg      = INTERFACETEST_ARG,
    udot     = INTERFACETEST_UDOT,
    rp       = INTERFACETEST_RP,
    e        = INTERFACETEST_E,
    alpha    = INTERFACETEST_ALPHA,
    fbcv     = INTERFACETEST_FBCV,
    sampling = INTERFACETEST_SAMPLING,
    deltaT   = INTERFACETEST_DELTAT, 
    length   = INTERFACETEST_LENGTH;




  UINT4 
    order  = INTERFACETEST_ORDER, 
    start=0,
    tail=0;
  

  UINT4 i,
    nspin=0;
    


  REAL8
    fdata[128];

  params->method = 0;



  for(i = 0; i < argc; i++)
    {
      if (strcmp(argv[i],"-approximant")==0)
	{
	  i++;
	  if (strcmp(argv[i],"PPN")==0)
	    params->method = PPN;
	  else if (strcmp(argv[i],"EOB")==0)
	      params->inspiral.approximant = EOB;
	  else if (strcmp(argv[i],"EOB")==0)
	    params->inspiral.approximant = EOB;
	  else if ((strcmp(argv[i],"TaylorCW")==0)  
	      || (strcmp(argv[i],"TAYLORCW")==0) 
	      || (strcmp(argv[i],"TCW")==0) 
	      )
	    params->method = TaylorCW;
	  else if ((strcmp(argv[i],"SPINORBITCW")==0)  
	      || (strcmp(argv[i],"SpinOrbitCW")==0) 
	      || (strcmp(argv[i],"SOCW")==0) 
	      )
	    params->method = SpinOrbitCW;
	  else if (strcmp(argv[i], "TaylorT1")==0)
	    params->inspiral.approximant = TaylorT1;
	  else if (strcmp(argv[i], "TaylorT2")==0)
	    params->inspiral.approximant = TaylorT2;
	  else if (strcmp(argv[i], "TaylorT3")==0)
	    params->inspiral.approximant = TaylorT3;
	  /*	  else if (strcmp(argv[i], "BCV")==0)
		  params->inspiral.approximant = BCV;*/
	  else
	    {
	      printf("Approximant not implemented\n");
	      break;
	    }




	}
      else if(strcmp(argv[i],"-tc")==0)
	tc = atof(argv[++i]);
      else if(strcmp(argv[i],"-m1")==0)
	m1 = atof(argv[++i]);
      else if(strcmp(argv[i],"-m2")==0)
	m2 = atof(argv[++i]);
      else if(strcmp(argv[i],"-order")==0)
	order = atoi(argv[++i]);
      else if(strcmp(argv[i],"-d")==0)
	d = atof(argv[++i]);
      else if(strcmp(argv[i],"-inclination")==0)
	inc = atof(argv[++i]);
      else if(strcmp(argv[i],"-position")==0)	
	{
	  ra  = atof(argv[++i]);
	  dec = atof(argv[++i]);
	}
      else if(strcmp(argv[i],"-psi")==0)
	psi = atof(argv[++i]);
      else if(strcmp(argv[i],"-phic")==0)
	phic = atof(argv[++i]);
      else if(strcmp(argv[i],"-freqlim")==0)	
	{
	  fi  = atof(argv[++i]);
	  ff = atof(argv[++i]);
	}
      else if(strcmp(argv[i],"-zeta2")==0)
	zeta2 = atof(argv[++i]);
      else if(strcmp(argv[i],"-omegas")==0)
	omegas = atof(argv[++i]);
      else if(strcmp(argv[i],"-theta")==0)
	theta = atof(argv[++i]);
      else if(strcmp(argv[i],"-t0")==0)
	t0 = atof(argv[++i]);
      else if(strcmp(argv[i],"-tcw")==0)	
	{
	  t1  = atof(argv[++i]);
	  t2  = atof(argv[++i]);
	}
      else if(strcmp(argv[i],"-sampling")==0)
	sampling = atof(argv[++i]);
      else if(strcmp(argv[i],"-deltat")==0)
	deltaT = atof(argv[++i]);
      else if(strcmp(argv[i],"-a")==0)
	{
	  a1  = atof(argv[++i]);
	  a2  = atof(argv[++i]);
	}
      else if(strcmp(argv[i],"-phi0")==0)
	phi0 = atof(argv[++i]);
      else if(strcmp(argv[i],"-f0")==0)
	{
	  f0 = atof(argv[++i]);
	  i++;
	  while ((*argv[i]!=45) && i <argc) 
	    {	      
	      fdata[nspin] = atof(argv[i]);
	      nspin++;
	      i++;
	    }	 
	  i--; 
	}      
      else if(strcmp(argv[i],"-arg")==0)
       arg = atof(argv[++i]);
      else if(strcmp(argv[i],"-udot")==0)
	udot = atof(argv[++i]);
      else if(strcmp(argv[i],"-rp")==0)
	rp = atof(argv[++i]);
      else if(strcmp(argv[i],"-e")==0)
	e = atof(argv[++i]);
      /* BCV */
      else if(strcmp(argv[i],"-alpha")==0)
	alpha = atof(argv[++i]);
      else if(strcmp(argv[i],"-fbcv")==0)
	fbcv = atof(argv[++i]);
      else if(strcmp(argv[i],"-start")==0)
	start = atoi(argv[++i]);
      else if(strcmp(argv[i],"-tail")==0)
	tail = atoi(argv[++i]);
      else if(strcmp(argv[i],"-length")==0)
	length = atoi(argv[++i]);
      else if(strcmp(argv[i],"-h")==0)
	{
	  /* Help message here*/
	  printf("help\n");
	  exit(0);
	}
      else
	{
	  fprintf(stderr,"%s: unrecognized option %s\n",
		  argv[0], argv[i]);

	  return 1;
	}
    }


 


  switch (params->inspiral.approximant)
    {
    case PPN:
      params->ppn.epoch.gpsSeconds   = params->ppn.epoch.gpsNanoSeconds = 0;
      params->ppn.position.latitude  = dec* LAL_PI/180.0;
      params->ppn.position.longitude = ra* LAL_PI/180.0;
      params->ppn.position.system    = COORDINATESYSTEM_EQUATORIAL;
      params->ppn.psi = psi* LAL_PI/180.0;

      params->ppn.mTot = m1 +m2;
      params->ppn.eta =   m1*m2   /( params->ppn.mTot*params->ppn.mTot );
      params->ppn.d = d * LAL_PC_SI*1e6;
      params->ppn.inc = inc * LAL_PI/180.0;
      params->ppn.phi = phic * LAL_PI/180.0;
      if (deltaT!=0) 
	params->ppn.deltaT = deltaT      ;
	  else
	params->ppn.deltaT = 1./ sampling;
      params->ppn.fStartIn = fi ;
      params->ppn.fStopIn = ff;;
      params->ppn.lengthIn = INTERFACETEST_LENGTHIN;  /* taille maximale du signal en points*/
      params->ppn.ppn = NULL;
      break;
    case SpinOrbitCW: 
      params->socw.epoch.gpsSeconds = params->socw.epoch.gpsNanoSeconds = 0;
      params->socw.spinEpoch.gpsSeconds = params->socw.spinEpoch.gpsNanoSeconds = 0;
      params->socw.orbitEpoch.gpsSeconds = params->socw.orbitEpoch.gpsNanoSeconds = 0;
      params->socw.position.latitude = dec* LAL_PI/180.0;
      params->socw.position.longitude = ra* LAL_PI/180.0;
      params->socw.position.system = COORDINATESYSTEM_EQUATORIAL;
      params->socw.psi = psi* LAL_PI/180.0;
      if (deltaT!=0) 
	params->socw.deltaT = deltaT;
      else
	params->socw.deltaT = 1./ sampling;
      params->socw.length = length;
      params->socw.aPlus=a1;
      params->socw.aCross=a2;
      params->socw.phi0=phi0 * LAL_PI/180.0;
      params->socw.f0=f0 ;
      params->socw.omega = arg * LAL_PI/180.0;
      params->socw.rPeriNorm = rp;   /* if =0 -> equivaut a Taylor*/
      params->socw.oneMinusEcc = 1 - e;
      params->socw.angularSpeed = udot ;  
      params->socw.f = NULL;
      
      if ( nspin ) {
	LALDCreateVector( status, &(params->socw.f), nspin );
	ClearStatus(status);
      }
      for (i=0; i<nspin; i++)	
	params->socw.f->data[i] = fdata[i];	  
	
      
      /* vp = rp * udot entre ]0, 1[
	 argument  = omega*/
      
      break;
    case TaylorCW:
      params->taylorcw.epoch.gpsSeconds = params->taylorcw.epoch.gpsNanoSeconds = 0;
      params->taylorcw.epoch.gpsSeconds = params->taylorcw.epoch.gpsNanoSeconds = 0;
      params->taylorcw.position.latitude = dec* LAL_PI/180.0;
      params->taylorcw.position.longitude = ra* LAL_PI/180.0;
      params->taylorcw.position.system = COORDINATESYSTEM_EQUATORIAL;
      params->taylorcw.psi = psi* LAL_PI/180.0;
      if (deltaT!=0) 
	params->taylorcw.deltaT = deltaT;
      else
	params->taylorcw.deltaT = 1./ sampling;
      params->taylorcw.length = length;
      params->taylorcw.aPlus=a1;
      params->taylorcw.aCross=a2;
      params->taylorcw.phi0=phi0 * LAL_PI/180.0;
      params->taylorcw.f0=f0 ;   /* ?? pouirquoi ce terme pi */
          
      params->taylorcw.f = NULL;
      if ( nspin  ) {
	LALDCreateVector( status,&(params->taylorcw.f), nspin );
	ClearStatus(status);
      }
      for (i=0; i<params->taylorcw.f->length; i++)
	params->taylorcw.f->data[i] = fdata[i];	  
      break;
    case TaylorT1:  case TaylorT2: case TaylorT3:
      params->inspiral.ieta=1;
      params->inspiral.mass1=m1;
      params->inspiral.mass2=m2;
      params->inspiral.startTime=0.0;
      params->inspiral.startPhase=0.0;
      params->inspiral.nStartPad=0;
      params->inspiral.nEndPad=0;
      params->inspiral.signalAmplitude = 1.; 
      params->inspiral.distance = d*LAL_PC_SI * 1e6; 
      params->inspiral.fLower=fi;
      params->inspiral.fCutoff=ff;
      params->inspiral.tSampling=sampling;
      params->inspiral.order=order;
      params->inspiral.inclination = inc *LAL_PI /180;      
      params->inspiral.massChoice=m1Andm2;
      if (deltaT!=0) 
	params->inspiral.tSampling = 1./ deltaT;
      else
	params->inspiral.tSampling = sampling;
      break;     
    case TaylorF1:
    case TaylorF2:
    case PadeF1:
    case BCV:
    case PadeT1:
      params->inspiral.ieta=1;
      params->inspiral.mass1=m1;
      params->inspiral.mass2=m2;
      params->inspiral.startTime=0.0;
      params->inspiral.startPhase=0.0;
      params->inspiral.nStartPad=0;
      params->inspiral.nEndPad=0;
      params->inspiral.signalAmplitude = 1.; 
      params->inspiral.distance = d*LAL_PC_SI * 1e6; 
      params->inspiral.fLower=fi;
      params->inspiral.fCutoff=ff;
      params->inspiral.tSampling=sampling;
      params->inspiral.order=order;
      params->inspiral.inclination = inc *LAL_PI /180;      
      params->inspiral.massChoice=m1Andm2;
      if (deltaT!=0) 
	params->inspiral.tSampling = 1./ deltaT;
      else
	params->inspiral.tSampling = sampling;
      params->inspiral.fendBCV = fbcv;
      params->inspiral.alpha   = alpha;
      break;
    case BCVSpin:
    case SpinTaylorT3:
    case EOB:
      params->inspiral.ieta=1;
      params->inspiral.mass1=m1;
      params->inspiral.mass2=m2;
      params->inspiral.startTime=0.0;
      params->inspiral.startPhase=0.0;
      params->inspiral.nStartPad=0;
      params->inspiral.nEndPad=0;
      params->inspiral.signalAmplitude = 1.; 
      params->inspiral.distance = d*LAL_PC_SI * 1e6; 
      params->inspiral.fLower=fi;
      params->inspiral.fCutoff=ff;
      params->inspiral.tSampling=sampling;
      params->inspiral.order=order;
      params->inspiral.inclination = inc *LAL_PI /180;      
      params->inspiral.massChoice=m1Andm2;
      if (deltaT!=0) 
	params->inspiral.tSampling = 1./ deltaT;
      else
	params->inspiral.tSampling = sampling;
      params->inspiral.OmegaS = omegas;
      params->inspiral.Theta  = theta;
      params->inspiral.Zeta2  = zeta2;
      break;
    default:      
      break;
  }




  return 0;
}






void
ClearStatus (
	     LALStatus   *status )
{
  if (status->statusPtr)
  {
    ClearStatus      (status->statusPtr);
    ClearStatus (status);
  }
}
