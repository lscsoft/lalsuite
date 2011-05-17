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
#include <lal/Window.h>
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

/**
\author Dietz, A. & Veitch, J.
\file
\ingroup inspiral
\brief Header file for the MCMC tools code.

\heading{Synopsis}
\code
#include <lal/LALInspiralMCMC.h>
\endcode

This header file covers routines that are used for the Markov Chain Monte Carlo algorithm tools.

\heading{Structures}

<ol>

<li> \c LALMCMCParameter:
Main structure that holds the different parameters that are used within the MCMC. The number, names or ranges is freely choosable for each of the parameter.

<dl>
<dt>\c tagLALMCMCParam</dt><dd> Pointer to a \c LALMCMCParam structure, which is a linked list over all parameters</dd>
<dt><tt>UINT4 dimension</tt></dt><dd> Dimension of the parameter space</dd>
<dt><tt>REAL8 logLikelihood</tt></dt><dd> The logarithm of the likelihood associated with this set of parameters </dd>
<dt><tt>REAL4 logPrior</tt></dt><dd> The logarithm of the prior associated with this set of parameters </dd>
</dl></li>

<li> \c LALMCMCParam:
Structure that contain the linked list structure and holds the current value.

<dl>
<dt><tt>LALMCMCParam next</tt></dt><dd> Pointer to the next \c LALMCMCParam structure (or a NULL pointer for the last element in this linked list).</dd>
<dt><tt>LALMCMCSubParam core</tt></dt><dd> Pointer to a \c LALMCMCSubParam structure that holds fixed values related to this parameter (see next structure).</dd>
<dt><tt>REAL8 value</tt></dt><dd> Actual value of this parameter</dd>
</dl></li>


<li> \c LALMCMCSubParam:
Structure that holds fixed properties for a single parameter.

<dl>
<dt><tt>char name</tt></dt><dd> Name of this parameter.</dd>
<dt><tt>REAL8 minVal</tt></dt><dd> Minimal allowed value for this parameter.</dd>
<dt><tt>REAL8 maxVal</tt></dt><dd> Maximal allowed value for this parameter.</dd>
<dt><tt>INT4 wrapping</tt></dt><dd> If set to 1, the value is being wrapped between \c minVal and \c maxVal (e.g. for any angle like a phase or right ascension).</dd>
<dt><tt>REAL4VECTOR chain</tt></dt><dd>  A \c REAL4Vector structure that holds the values of the chain for this parameter.</dd>
</dl></li>

<li> \c LALMCMCInput:
Structure that holds all data needed for the MCMC algorithm. It contains input data, output data, as well as flags and parameters used for the MCMC algorithm.

<dl>
<dt><tt>RandomParams randParams</tt></dt><dd> Parameter for random number generation.</dd>
<dt><tt>FindChirpFilterInput fcFilterInput</tt></dt><dd> A FindChirpFilterInput structure containing the input data</dd>
<dt><tt>FindChirpFilterParams fcFilterParams</tt></dt><dd> A FindChirpFilterParams structure containing the parameters</dd>
<dt><tt>FindChirpDataParams fcDataParams</tt></dt><dd> Parameters describing the data </dd>
<dt><tt>SnglInspiralTable     *inspiralTable</tt></dt><dd> A pointer to a single\_inspiral table containing parameters of the trigger</dd>
<dt><tt>MCMCInitFunction *funcInit</tt></dt><dd> A pointer to a function that initializes the parameter structure</dd>
<dt><tt>MCMCLikelihoodFunction *funcLikelihood</tt></dt><dd> A pointer to a function that calculates the logarithm of the likelihood</dd>
<dt><tt>MCMCPriorFunction *funcPrior</tt></dt><dd> A pointer to a function that calculates the logarithm of the prior</dd>
<dt><tt>InspiralTemplate tmpltPtr</tt></dt><dd> A InspiralTemplate structure to hold the parameters to create a template </dd>
<dt><tt>FindChirpTmpltParams fcTmpltParams</tt></dt><dd> A FindChirpTmpltParams structure to hold parameters for creating the template</dd>
<dt><tt>Approximant approximant</tt></dt><dd> The approximant used to filter the data (NA)</dd>
<dt><tt>UINT4 verbose</tt></dt><dd> The verbosity flag (NA)</dd>


<dt><tt>UINT4 dim</tt></dt><dd> Number of dimensions of the parameter space</dd>
<dt><tt>UINT4 counter</tt></dt><dd> Counter of the MCMC iteration in general</dd>
<dt><tt>UINT4 counterState</tt></dt><dd> Counter of the MCMC iteration within the current state</dd>
<dt><tt>UINT4 counterAccept</tt></dt><dd> Counter of the MCMC iteration within the accepted steps</dd>
<dt><tt>UINT4 counterAcceptDraw</tt></dt><dd> Counter of the MCMC iteration which are the drawn ones</dd>

<dt><tt>UINT4 numberDraw</tt></dt><dd> Number of values to be drawn after all the possible methods below came into place\\</dd>

<dt><tt>UINT4 useAnnealing</tt></dt><dd> Flag to activate annealing</dd>
<dt><tt>UINT4 numberAnneal</tt></dt><dd> Number of iterations used for matrix annealing.</dd>
<dt><tt>REAL4 annealTempBegin</tt></dt><dd> Initial annealing temperature</dd>
<dt><tt>REAL4 annealTemp</tt></dt><dd> Actual annealing temperature</dd>
<dt><tt>UINT4 annealingSteps</tt></dt><dd> TO BE SPECIFIED</dd>

<dt><tt>UINT4 useScaling</tt></dt><dd> Flag to activate the scaling method</dd>
<dt><tt>REAL4 scaling</tt></dt><dd> Initial scaling value (should be named Begin or so,. see annealing)</dd>
<dt><tt>REAL4 scalePeak</tt></dt><dd> Initial scaling value (e.g. 50.0)</dd>
<dt><tt>REAL4 scaleNormal</tt></dt><dd> Normal scaling value (1.0)</dd>
<dt><tt>REAL4 scaleQ</tt></dt><dd> an internal parameter </dd>
<dt><tt>REAL4 scalePA</tt></dt><dd> an internal parameter </dd>

<dt><tt>UINT4 flagAcceptRatio</tt></dt><dd> Flag to activate the acceptance-ratio method</dd>
<dt><tt>UINT4 acceptRatioCounter</tt></dt><dd> Factor for acceptance ratio method</dd>
<dt><tt>UINT4 acceptRatioNorm</tt></dt><dd> Norming factor for acceptance ratio method</dd>
<dt><tt>REAL4 acceptThreshold</tt></dt><dd> Threshold </dd>
<dt><tt>Approximant approximant</tt></dt><dd> Approximant to use</dd>

<dt><tt>UINT4 useUpdate</tt></dt><dd> Flag to activate matrix updating.</dd>
<dt><tt>UINT4 updateNumber</tt></dt><dd> Number of iterations used to update the covariance matrix.</dd>
<dt><tt>UINT4 updateOffset</tt></dt><dd> Offset value used in matrix updating</dd>
<dt><tt>UINT4 updateCounter</tt></dt><dd> Internal counter used by the updating algorithm</dd>
<dt><tt>REAL8* mean</tt></dt><dd> A vector containing the mean values of each parameter</dd>
<dt><tt>REAL8* xdiff</tt></dt><dd> A vector used for updating the matrix</dd>
<dt><tt>REAL8* ndiff</tt></dt><dd> A vector used for updating the matrix</dd>


<dt><tt>UINT4 flagBurnin</tt></dt><dd> Flag to activate the burn-in method</dd>
<dt><tt>UINT4 burninNumber</tt></dt><dd> Minimum steps after which the burn-in period is checked</dd>
<dt><tt>UINT4 burninStep</tt></dt><dd> The number of steps between two checks for burn-in</dd>
<dt><tt>UINT4 burninCounter</tt></dt><dd> Internal counter</dd>
<dt><tt>UINT4 burninTime</tt></dt><dd> The number of iteration will be stored after the burn-in is reached</dd>
<dt><tt>REAL4 burninThreshold</tt></dt><dd> Threshold for this chain of having reached burnin.

</dd>
</dl>

</li>
</ol>

*/


NRCSID( LALINSPIRALMCMCH, "$Id: LALInspiralMCMC.h,v 1.79 2007/02/19 15:52:17 thomas Exp $" );

/**\name Error Codes */ /*@{*/
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
/*@}*/



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




typedef struct
tagLALMCMCSubParam
{
  char        name[30];
  REAL8       minVal;
  REAL8       maxVal;
  INT4        wrapping; /* 0=no, 1=yes, -1=fixed */
  REAL4Vector *chain;

}  LALMCMCSubParam;




typedef struct
tagLALMCMCParam
{
  struct tagLALMCMCParam    *next;
  struct tagLALMCMCSubParam *core;
  REAL8                     value;
}  LALMCMCParam;




typedef struct
tagLALMCMCParameter
{
  struct tagLALMCMCParam* param;
  UINT4                   dimension;
  REAL8                   logLikelihood;
  REAL8                   logPrior;
}  LALMCMCParameter;



typedef void (MCMCInitFunction)(
  LALMCMCParameter  *parameter,
  void *input);


struct tagLALMCMCInput;


typedef REAL8 (MCMCLikelihoodFunction)(
    struct tagLALMCMCInput *inputMCMC,
    LALMCMCParameter  *parameter);



typedef REAL8 (MCMCPriorFunction)(
   struct tagLALMCMCInput      *inputMCMC,
   LALMCMCParameter  *parameter);



typedef struct
tagLALMCMCInput
{
  RandomParams *randParams;

  UINT4                     numberDataStreams;
  CHAR*                     ifoID[MAXDET];
  CHAR*						dumpfile; /* Likelihod function should dump data if this is not null */
  REAL8TimeSeries*          segment[MAXDET];
  REAL8FrequencySeries*     invspec[MAXDET];
  COMPLEX16FrequencySeries* stilde[MAXDET];
  LALDetector*				detector[MAXDET];

  SnglInspiralTable         *inspiralTable;
  SimInspiralTable			*injectionTable;
  REAL8FFTPlan *fwdplan;
  REAL8FFTPlan *revplan;
  REAL4FFTPlan *likelihoodPlan;
  REAL4FFTPlan *likelihoodRevPlan;
  REAL8Window *window; /* Window for FFTing the data */
  LIGOTimeGPS epoch;
  REAL4   fLow;
  REAL8   deltaT;
  REAL8   deltaF;
  UINT4   numseg;
  UINT4   numPoints; /* numPoints */
  UINT4   stride;   /* ovrlap */
  Approximant approximant; /* Approximant to use for this model */
  INT4	  ampOrder; /* Amplitude order to use with Higher Harmonic waveforms */
	                /* Setting = 0 means Newtonian amplitude */
  LALPNOrder phaseOrder;
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

/* For plus and cross polarisations in PhenSpinRD */
  REAL4Vector* Fwfp;
  REAL4Vector* Fwfc;
  REAL4FFTPlan *longplan;
  UINT4 mylength;
}  LALMCMCInput;





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

INT4 XLALMCMCJumpHarmonic(
  LALMCMCInput *inputMCMC,
  LALMCMCParameter *parameter
     );

void XLALMCMCJumpSingle(
  LALMCMCInput *inputMCMC,
  LALMCMCParameter *parameter,
  gsl_matrix       *covMat
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

INT4 XLALMCMCCheckWrapping(LALMCMCParameter *parameter,
						   const char *name);

int PriorIsSane(LALMCMCParameter *parameter);

#ifdef  __cplusplus
}
#endif

#endif /* _LALINSPIRAL_H */
