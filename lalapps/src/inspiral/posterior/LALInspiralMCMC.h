/* <lalVerbatim file="LALInspiralMCMCHV">
Author: Dietz, A. & Veitch, J.
$Id: LALInspiralMCMC.h,v 1.79 2007/02/19 15:52:17 thomas Exp $
</lalVerbatim>  */


/*  <lalLaTeX>

\section{Header \texttt{LALInspiralMCMC.h}}
\label{s:LALInspiralMCMC.h}

Header file for the MCMC tools code.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALInspiralMCMC.h>
\end{verbatim}

\noindent This header file covers routines that are used for the Markov Chain Monte Carlo algorithm tools.

</lalLaTeX> */

#ifndef _LALINSPIRALMCMC_H
#define _LALINSPIRALMCMC_H

# include <math.h>
# include <stdio.h>
# include <stdlib.h>


# include <lal/LALStdlib.h>
# include <lal/LALConstants.h>
# include <lal/SimulateCoherentGW.h>
# include <lal/GeneratePPNInspiral.h>
# include <lal/LIGOMetadataTables.h>
# include <lal/LALDatatypes.h>
# include <lal/FindChirp.h>
#include <lal/LALDetectors.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>


#ifdef  __cplusplus
extern "C" { 
#endif

NRCSID( LALINSPIRALMCMCH, "$Id: LALInspiralMCMC.h,v 1.79 2007/02/19 15:52:17 thomas Exp $" );

/*  <lalLaTeX>
\subsection*{Error codes}
</lalLaTeX>  */

/* <lalErrTable> */
#define LALINSPIRALH_ENULL           1
#define LALINSPIRALH_EMEM            2
#define LALINSPIRALH_EDIV0           3
#define LALINSPIRALH_ESIZE           4
#define LALINSPIRALH_ECHOICE         5
#define LALINSPIRALH_EORDER          6 
#define LALINSPIRALH_EAPPROXIMANT    7 

#define LALINSPIRALH_MSGENULL         "Arguments contained an unexpected null pointer"
#define LALINSPIRALH_MSGEMEM          "Memory allocation error"
#define LALINSPIRALH_MSGEDIV0         "Division by zero"
#define LALINSPIRALH_MSGESIZE         "Invalid input range"
#define LALINSPIRALH_MSGECHOICE       "Invalid choice for an input parameter"
#define LALINSPIRALH_MSGEORDER        "unknown order specified"
#define LALINSPIRALH_MSGEAPPROXIMANT  "Invalid model"
#define MAXDET 5			/* Maximum number of data streams/detectors to accept */
/** ---------------------------------------------------------------------  
</lalErrTable> */



/* <lalLaTeX>

\section*{Structures}

\begin{enumerate}

\item \texttt{LALMCMCParameter:}
Main structure that holds the different parameters that are used within the MCMC. The number, names or ranges is freely choosable for each of the parameter.
\input{LALMCMCParameterH}
\begin{description}
\item[\texttt{tagLALMCMCParam}] Pointer to a {\tt LALMCMCParam} structure, which is a linked list over all parameters
\item[\texttt{UINT4 dimension}] Dimension of the parameter space
\item[\texttt{REAL8 logLikelihood}] The logarithm of the likelihood associated with this set of parameters 
\item[\texttt{REAL4 logPrior}] The logarithm of the prior associated with this set of parameters 
\end{description}

\item \texttt{LALMCMCParam:}
Structure that contain the linked list structure and holds the current value. 
\input{LALMCMCParamH}
\begin{description}
\item[\texttt{LALMCMCParam next}] Pointer to the next {\tt LALMCMCParam} structure (or a NULL pointer for the last element in this linked list).
\item[\texttt{LALMCMCSubParam core}] Pointer to a {\tt LALMCMCSubParam} structure that holds fixed values related to this parameter (see next structure).
\item[\texttt{REAL8 value}] Actual value of this parameter
\end{description}


\item \texttt{LALMCMCSubParam:}
Structure that holds fixed properties for a single parameter.
\input{LALMCMCSubParamH}
\begin{description}
\item[\texttt{char name}] Name of this parameter.
\item[\texttt{REAL8 minVal}] Minimal allowed value for this parameter.
\item[\texttt{REAL8 maxVal}] Maximal allowed value for this parameter.
\item[\texttt{INT4 wrapping}] If set to 1, the value is being wrapped between {\tt minVal} and {\tt maxVal} (e.g. for any angle like a phase or right ascension).
\item[\texttt{REAL4VECTOR chain}]  A {\tt REAL4Vector} structure that holds the values of the chain for this parameter.
\end{description}

\item \texttt{LALMCMCInput:}
Structure that holds all data needed for the MCMC algorithm. It contains input data, output data, as well as flags and parameters used for the MCMC algorithm.
\input{LALMCMCInputH}
\begin{description}
\item[\texttt{RandomParams randParams}] Parameter for random number generation.
\item[\texttt{FindChirpFilterInput fcFilterInput}] A FindChirpFilterInput structure containing the input data
\item[\texttt{FindChirpFilterParams fcFilterParams}] A FindChirpFilterParams structure containing the parameters
\item[\texttt{FindChirpDataParams fcDataParams}] Parameters describing the data 
\item[\texttt{SnglInspiralTable     *inspiralTable}] A pointer to a single\_inspiral table containing parameters of the trigger
\item[\texttt{MCMCInitFunction *funcInit}] A pointer to a function that initializes the parameter structure
\item[\texttt{MCMCLikelihoodFunction *funcLikelihood}] A pointer to a function that calculates the logarithm of the likelihood
\item[\texttt{MCMCPriorFunction *funcPrior}] A pointer to a function that calculates the logarithm of the prior
\item[\texttt{InspiralTemplate tmpltPtr}] A InspiralTemplate structure to hold the parameters to create a template 
\item[\texttt{FindChirpTmpltParams fcTmpltParams}] A FindChirpTmpltParams structure to hold parameters for creating the template
\item[\texttt{Approximant approximant}] The approximant used to filter the data (NA)
\item[\texttt{UINT4 verbose}] The verbosity flag (NA)


\item[\texttt{UINT4 dim}] Number of dimensions of the parameter space
\item[\texttt{UINT4 counter}] Counter of the MCMC iteration in general
\item[\texttt{UINT4 counterState}] Counter of the MCMC iteration within the current state
\item[\texttt{UINT4 counterAccept}] Counter of the MCMC iteration within the accepted steps
\item[\texttt{UINT4 counterAcceptDraw}] Counter of the MCMC iteration which are the drawn ones

\item[\texttt{UINT4 numberDraw}] Number of values to be drawn after all the possible methods below came into place\\

\item[\texttt{UINT4 useAnnealing}] Flag to activate annealing
\item[\texttt{UINT4 numberAnneal}] Number of iterations used for matrix annealing.
\item[\texttt{REAL4 annealTempBegin}] Initial annealing temperature
\item[\texttt{REAL4 annealTemp}] Actual annealing temperature
\item[\texttt{UINT4 annealingSteps}] TO BE SPECIFIED

\item[\texttt{UINT4 useScaling}] Flag to activate the scaling method
\item[\texttt{REAL4 scaling}] Initial scaling value (should be named Begin or so,. see annealing)
\item[\texttt{REAL4 scalePeak}] Initial scaling value (e.g. 50.0)
\item[\texttt{REAL4 scaleNormal}] Normal scaling value (1.0)
\item[\texttt{REAL4 scaleQ}] an internal parameter 
\item[\texttt{REAL4 scalePA}] an internal parameter 

\item[\texttt{UINT4 flagAcceptRatio}] Flag to activate the acceptance-ratio method
\item[\texttt{UINT4 acceptRatioCounter}] Factor for acceptance ratio method
\item[\texttt{UINT4 acceptRatioNorm}] Norming factor for acceptance ratio method
\item[\texttt{REAL4 acceptThreshold}] Threshold 
\item[\texttt{Approximant approximant}] Approximant to use

\item[\texttt{UINT4 useUpdate}] Flag to activate matrix updating.
\item[\texttt{UINT4 updateNumber}] Number of iterations used to update the covariance matrix.
\item[\texttt{UINT4 updateOffset}] Offset value used in matrix updating
\item[\texttt{UINT4 updateCounter}] Internal counter used by the updating algorithm
\item[\texttt{REAL8* mean}] A vector containing the mean values of each parameter
\item[\texttt{REAL8* xdiff}] A vector used for updating the matrix
\item[\texttt{REAL8* ndiff}] A vector used for updating the matrix


\item[\texttt{UINT4 flagBurnin}] Flag to activate the burn-in method
\item[\texttt{UINT4 burninNumber}] Minimum steps after which the burn-in period is checked
\item[\texttt{UINT4 burninStep}] The number of steps between two checks for burn-in
\item[\texttt{UINT4 burninCounter}] Internal counter
\item[\texttt{UINT4 burninTime}] The number of iteration will be stored after the burn-in is reached
\item[\texttt{REAL4 burninThreshold}] Threshold for this chain of having reached burnin.


\end{description}


\end{enumerate}

--------------------------------------------------------------------- </lalLaTeX>  */

/* 
  prototypes for MCMC 
*/

/** enum containing the different ways in which the mcmc is set up **/
typedef enum
{
  unknownMode,
  modeTest,
  modeEOB,
  modeSpinning, 
  modeTaylor
}
MCMCmode;



/* <lalVerbatim file="LALMCMCSubParamH">  */
typedef struct
tagLALMCMCSubParam
{
  char        name[30]; 
  REAL8       minVal;  
  REAL8       maxVal;  
  INT4        wrapping;
  REAL4Vector *chain;

}  LALMCMCSubParam;
/* </lalVerbatim>  */


/* <lalVerbatim file="LALMCMCParamH">  */
typedef struct
tagLALMCMCParam
{
  struct tagLALMCMCParam    *next; 
  struct tagLALMCMCSubParam *core; 
  REAL8                     value;    
}  LALMCMCParam;
/* </lalVerbatim>  */


/* <lalVerbatim file="LALMCMCParameterH">  */
typedef struct
tagLALMCMCParameter
{
  struct tagLALMCMCParam* param;
  UINT4                   dimension;      
  REAL8                   logLikelihood; 
  REAL8                   logPrior;
}  LALMCMCParameter;
/* </lalVerbatim>  */

/* <lalVerbatim file="MCMCInitFunctionH">  */
typedef void (MCMCInitFunction)(
  LALMCMCParameter  *parameter,
  void *input);
/* </lalVerbatim>  */

struct tagLALMCMCInput;

/* <lalVerbatim file="MCMCLikelihoodFunctionH">  */
typedef REAL8 (MCMCLikelihoodFunction)(
    struct tagLALMCMCInput *inputMCMC,
    LALMCMCParameter  *parameter);
/* </lalVerbatim>  */

/* <lalVerbatim file="MCMCPriorFunctionH">  */
typedef REAL8 (MCMCPriorFunction)(
   struct tagLALMCMCInput      *inputMCMC,
   LALMCMCParameter  *parameter);
/* </lalVerbatim>  */

/* <lalVerbatim file="LALMCMCInputH">  */
typedef struct
tagLALMCMCInput
{ 
  RandomParams *randParams;
   
  UINT4                     numberDataStreams;
  CHAR*                     ifoID[MAXDET];
  REAL8TimeSeries*          segment[MAXDET];
  REAL8FrequencySeries*     invspec[MAXDET];
  COMPLEX16FrequencySeries* stilde[MAXDET];
  LALDetector*				detector[MAXDET];

  SnglInspiralTable         *inspiralTable;
  SimInspiralTable			*injectionTable;
  REAL8FFTPlan *fwdplan;
  REAL8FFTPlan *revplan;
  REAL8Window *window; /* Window for FFTing the data */
  LIGOTimeGPS epoch;
  REAL4   fLow;
  REAL8   deltaT;
  REAL8   deltaF;
  UINT4   numseg;
  UINT4   numPoints; /* numPoints */
  UINT4   stride;   /* ovrlap */
  Approximant approximant; /* Approximant to use for this model */

  MCMCmode mode;
  MCMCLikelihoodFunction *funcLikelihood; /* engine for likelihood */
  MCMCInitFunction       *funcInit;       /* engine for init function */
  MCMCPriorFunction      *funcPrior;      /* engine for prior */

  UINT4 verbose;    /* verbosity flag */

  UINT4 dim;           /* dimension of the parameter space */
  UINT4 counter;      /* overall counter of the actual iteration */
  UINT4 counterState; /* counter of the iteration within current state */
  UINT4 counterAccept; /* overal counter of accetped steps */
  UINT4 counterAcceptDraw; 

  /* some user arguments */
  UINT4 numberDraw;

  /* parameters related to annealing */
  UINT4 useAnnealing;   /* switch for usinh annealing */
  UINT4 numberAnneal;
  REAL4 annealingTempBegin;  /* starting annealing temperature */
  REAL4 annealingTemp ;       /* actual annealing temperature */
  UINT4 annealingSteps;      /* X */

  /* parameters related to scaling */
  UINT4 useScaling;
  REAL4 scaling;         /* actual scaling factor */
  REAL4 scalePeak;       /* peaking scaling factor */
  REAL4 scaleNormal;     /* normal scaling factor, should be 1.0 */
  REAL4 scaleQ;        
  REAL4 scalePA;

  /* parameters related to covariance updating */
  UINT4 useUpdate;         /* switch for using matrix updating */
  UINT4 updateNumber ;     /* number of draws for matrix updating */
  UINT4 updateOffset ;
  UINT4 updateCounter ;    /* internal counter used for the updating process */
  REAL8* mean;
  REAL8* xdiff;
  REAL8* ndiff;

  /* parameters related to burn-in */
  UINT4 flagBurnin;      /* switch to check if burn-in period is reached */
  UINT4 burninNumber;   /* minimum steps for burnin checking */
  UINT4 burninStep;     /* The number of steps used to check burn-in */
  UINT4 burninCounter;   /* Internal counter used */
  UINT4 burninTime ;     /* step at which the chain is 'burnt-in */
  REAL4 burninThreshold;
  UINT4 burninMaxNumber;  /* maximum number of trials */

	/* Parameter for nested sampling */
	UINT4 Nlive;
	LALMCMCParameter **Live;

}  LALMCMCInput;
/* </lalVerbatim>  */




typedef enum
{
  unknownState,
  doScaling,
  doAnnealing, 
  doUpdating,
  doDrawing
}
StateMCMC;





/* Function prototypes */

/* --- MCMC code ---- */


/*  <lalLaTeX>
\newpage\input{LALInspiralMCMCC}
</lalLaTeX>  */


void 
printMatrix( gsl_matrix *covMat, int dim);

void 
printState( StateMCMC mode );




void 
XLALMCMCBasicMetro(
   LALMCMCParameter **parameter,
   LALMCMCInput *inputMCMC);

UINT4 
XLALMCMCBasicSample(
   LALMCMCInput *inputMCMC,
   LALMCMCParameter **pparameter, 
   REAL4 *posterior);

void 
XLALMCMCBasicJump(
   LALMCMCInput *inputMCMC,
   LALMCMCParameter *parameter);


void 
XLALMCMCMetro(
   LALMCMCParameter **parameter,
   LALMCMCInput *inputMCMC);

void 
XLALMCMCCheckAnnealing(
  gsl_matrix *covMat,
  gsl_matrix *inputMat,
  LALMCMCInput *inputMCMC);  

INT4
XLALMCMCCheckBurnin(
  LALMCMCInput *inputMCMC,  
  LALMCMCParameter *parameter);

void
XLALMCMCCheckUpdate(
   LALMCMCInput *inputMCMC,
   LALMCMCParameter *parameter,
   gsl_matrix *covMat);

void
XLALMCMCCheckAcceptRatio(
  LALMCMCInput *inputMCMC,
  int move);

UINT4 
XLALMCMCSample(
   LALMCMCInput *inputMCMC,
   LALMCMCParameter **pparameter, 
   REAL4 *posterior, 
   gsl_matrix *covMat);

void 
XLALMCMCJump(
   LALMCMCInput *inputMCMC,
   LALMCMCParameter *parameter,
   gsl_matrix *covMat);

void 
XLALMCMCJumpIntrinsic(
  LALMCMCInput     *inputMCMC,
  LALMCMCParameter *parameter, 
  gsl_matrix       *covMat
  );

void XLALMCMCCyclicReflectiveBound(LALMCMCParameter *parameter);


void XLALMCMCGetCartesianPos(REAL8 vec[3],LALMCMCParameter *parameter);

void CartesianToSkyPos(REAL8 pos[3],LALMCMCParameter *parameter);


void crossProduct(REAL8 out[3],REAL8 x[3],REAL8 y[3]);

INT4 XLALMCMCDifferentialEvolution(
        LALMCMCInput *inputMCMC,
        LALMCMCParameter *parameter);

INT4 XLALMCMCReflectDetPlane(
	LALMCMCInput *inputMCMC,
	LALMCMCParameter *parameter);


void XLALMCMCRotateSky(
	LALMCMCInput *inputMCMC,
	LALMCMCParameter *parameter
	);

int XLALMCMC1PNMasseta(LALMCMCInput *inputMCMC, LALMCMCParameter *parameter);

INT4 XLALMCMCCheckParameter(
			   LALMCMCParameter *parameter,
			   const char *name);

void
XLALMCMCAddParam(
   LALMCMCParameter   *parameter,
   const char         *name, 
   REAL8              value,
   REAL8              minValue,
   REAL8              maxValue,
   INT4               wrapping);

LALMCMCParam*
XLALMCMCGetParam(
   LALMCMCParameter* parameter,
   const char *name);

REAL8
XLALMCMCGetParameter(
   LALMCMCParameter* parameter,
   const char *name);


void
XLALMCMCSetParameter(
   LALMCMCParameter* parameter,
   const char* name,
   REAL8 value);


void
XLALMCMCCopyPara(
   LALMCMCParameter **parameterOut,
   LALMCMCParameter *parameterIn);
       

void
XLALMCMCFreePara(
    LALMCMCParameter *parameter);


void
XLALMCMCDestroyPara(
    LALMCMCParameter **parameter);

void
XLALMultiStudentDeviates( 
   REAL4Vector  *vector,
   gsl_matrix   *matrix,
   UINT4         dim,
   UINT4         n,
   RandomParams *randParam);


void
XLALMultiNormalDeviates( 
   REAL4Vector  *vector,
   gsl_matrix   *matrix,
   UINT4         dim,
   RandomParams *randParam);


UINT4
XLALCheckPositiveDefinite( 
   gsl_matrix       *matrix,
   UINT4         dim);


#ifdef  __cplusplus
}
#endif

#endif /* _LALINSPIRAL_H */
