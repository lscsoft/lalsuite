/*  <lalVerbatim file="LALInspiralInitCV">
Author: Cokelaer T.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralInit.c}}
Module to initialize some parameters for waveform generation. 

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralRestrictedInitCP}
\idx{LALInspiralRestrictedInit()}

\subsubsection*{Description}
The input parameters is an InspiralTemplate structure which provides the waveform parameters
such as masses, lower frequency\dots The function \texttt{LALInspiralInit} calls the 
\texttt{LALInspiralParameterCalc} function in order to  compute all the mass parameters. Then, 
\texttt{LALInspiralRestrictedAmplitude} function is called to get the restricted newtonian 
amplitude. LALInspiralWavelength, LALInspiralSetup and LALInspiralChooseModel are also called 
in order to estimate the waveform length which is stored in an output structure called 
\texttt{InspiralInit}. We also stored Energy, flux and evolution function of flux and energy in 
that structure.

The  \texttt{LALInspiralChooseModel} function might failed or send a non zero status code. 
That function force it to be zero therefore the codes which  use LALInspiralInit (mainly 
injection code right now) won't stopped. Of course, if status code is non zero, we have to keep 
trace of it. Thus, the length of the waveform is fixed to zero in case of problems such as
negative length, cutoff frequency lower than the lower cutoff frequency \dots. 

\subsubsection*{Uses}
\texttt{LALInspiralParameterCalc}\\
\noindent\texttt{LALInspiralRestrictedAmplitude}\\
\noindent\texttt{LALInspiralWaveLength}
\noindent\texttt{LALInspiralChooseModel}
\noindent\texttt{LALInspiralSetup}

\subsubsection*{Notes}
There is only one assert on the InspiralTemplate variable since  all relevant asserts
are already included in the different functions which are called throughout the LALInspiralInit
function.
\vfill{\footnotesize\input{LALInspiralInitCV}}
</lalLaTeX>  */


#include <lal/LALInspiral.h>
#define  LALINSPIRALINIT_LENGTHOVERESTIMATION  0.1       /* 10 % */
#define  LALINSPIRALINIT_MINIMALWAVELENGTH     0.03125   /* 64 bins with a 2048Hz sampling*/

NRCSID (LALINSPIRALAMPLITUDEC, "$Id$");

/*  <lalVerbatim file="LALInspiralRestrictedInitCP"> */
void 
LALInspiralInit (LALStatus        *status, 
		 InspiralTemplate *params, 
		 InspiralInit     *paramsInit)
{ /* </lalVerbatim> */
  
  UINT4 ndx;
  REAL8 x;
  CHAR message[256];

  INITSTATUS (status, "LALInspiralInit", LALINSPIRALAMPLITUDEC );
  ATTATCHSTATUSPTR(status);

  ASSERT( params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );


  LALInspiralParameterCalc(status->statusPtr,  params);
  CHECKSTATUSPTR(status); 
  
  LALInspiralRestrictedAmplitude(status->statusPtr, params);
  CHECKSTATUSPTR(status);
  
  LALInspiralSetup(status->statusPtr, &(paramsInit->ak), params);
  CHECKSTATUSPTR(status);
  
  LALInspiralChooseModel(status->statusPtr, &(paramsInit->func), &(paramsInit->ak), params);

  /* The parameters have been initialized now. However, we can have some problems 
     with the LALInspiralChooseModel related to bad estimation of the length. 

     We first need to check that the length is not wrong. 
     Then, to check that flso is > fLow

     keep a security length higher than the one given by ChooseModel
  */

  if( params->fCutoff < params->fLower){
    LALWarning(status,  LALINSPIRALH_MSGEFLOWER);
    status->statusPtr->statusCode = 0;    
    paramsInit->nbins = 0;

    sprintf(message, "#Estimated Length (seconds) = %f | fCutoff = %f",
	    paramsInit->ak.tn, params->fCutoff);
    LALInfo(status, message);

    
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }  
  
  if( paramsInit->ak.tn <=0 || params->tC <= 0){
    LALWarning(status,  LALINSPIRALH_MSGESIZE);
    status->statusPtr->statusCode = 0;    
    paramsInit->nbins = 0;
    sprintf(message, "#Estimated Length (seconds) = %f ",
	    paramsInit->ak.tn);
    LALInfo(status, message);
    
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }  
  
  if (status->statusPtr->statusCode == 0){/* if everything is fine is ChooseModel then we
					     estimate the waveform length. */
    /*we add a minimal value and 10 % of overestimation */
    x	= (1.+ LALINSPIRALINIT_LENGTHOVERESTIMATION) 
      * (paramsInit->ak.tn + LALINSPIRALINIT_MINIMALWAVELENGTH ) * params->tSampling 
      + params->nStartPad + params->nEndPad;
    ndx 	= ceil(log10(x)/log10(2.));
    paramsInit->nbins =  pow(2, ndx);
 

    /*now we can free memory */
    CHECKSTATUSPTR(status);

    sprintf(message, "#Estimated Length (seconds) = %f | Allocated length (bins) = %d", 
	    paramsInit->ak.tn, 
	    paramsInit->nbins);
    LALInfo(status, message);

    
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }
  else { /*otherwise size is zero */
    sprintf(message,
	    "Can't get size of the following waveform: totalMass = %f, fLower = %f, approximant = %d @ %fPN"
	    , params->mass1 + params->mass2, params->fLower, params->approximant, params->order/2.);
    LALWarning(status, message);

    status->statusPtr->statusCode = 0;
    paramsInit->nbins = 0;

    /*now we can free memory */
    CHECKSTATUSPTR(status);
    
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }
  
}

