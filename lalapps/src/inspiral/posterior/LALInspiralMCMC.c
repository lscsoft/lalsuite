/*  <lalVerbatim file="LALInspiralMCMCCV">
Author: Sathyaprakash, B. S.
$Id: LALInspiralPhase.c,v 1.9 2003/04/14 00:27:22 sathya Exp $
</lalVerbatim>  */


#if 0
/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralMCMC.c} }

The file \texttt{LALInspiralMCMC} contains tools to perform a Monte Carlo Markov Chain parameter estimation computation on gravitational wave data. 

\subsubsection*{Prototypes}
\vspace{0.1in}


\input{XLALMCMCMetroCP}
\index{\verb&XLALMCMCMetro()&}

\input{XLALMCMCCheckAnnealingCP}
\index{\verb&XLALMCMCCheckAnnealing()&}

\input{XLALMCMCCheckBurninCP}
\index{\verb&XLALMCMCCheckBurnin()&}

\input{XLALMCMCCheckUpdateCP}
\index{\verb&XLALMCMCCheckUpdate()&}

\input{XLALMCMCCheckAcceptRatioCP}
\index{\verb&XLALMCMCCheckAcceptRatio()&}

\input{XLALMCMCSampleCP}
\index{\verb&XLALMCMCSample()&}

\input{XLALMCMCJumpCP}
\index{\verb&XLALMCMCJump()&}


\input{XLALMCMCAddParamCP}
\index{\verb&XLALMCMCAddParam()&}

\input{XLALMCMCGetParamCP}
\index{\verb&XLALMCMCGetParam()&}

\input{XLALMCMCGetParameterCP}
\index{\verb&XLALMCMCGetParameter()&}

\input{XLALMCMCSetParameterCP}
\index{\verb&XLALMCMCSetParameter()&}

\input{XLALMCMCCopyParaCP}
\index{\verb&XLALMCMCCopyPara()&}

\input{XLALMCMCFreeParaCP}
\index{\verb&XLALMCMCFreePara()&}

\input{XLALMCMCDestroyParaCP}
\index{\verb&XLALMCMCDestroyPara()&}

\input{XLALMultiStudentDeviatesCP}
\index{\verb&XLALMultiStudentDeviates()&}

\input{XLALMultiNormalDeviatesCP}
\index{\verb&XLALMultiNormalDeviates()&}

\input{XLALCheckPositiveDefiniteCP}
\index{\verb&XLALCheckPositiveDefinite()&}


\subsubsection*{Description}

This package contains routines needed for doing a MCMC calculation within the LAL framework. The only functions called from outside this package is to {\tt XLALMCMCMetro}. 

The function {\tt XLALMCMCMetro} is the main function to call to start the actual Markov Chain given the parameters in the \texttt{MCMCParameter} structure as a set of starting value. This function performs a Metropolis Hasting sampling. 

The function {\tt XLALMCMCSample} samples the next element of the chain, i.e. draws a new proposal parameter set using {\tt XLALMCMCJump} and apply the acceptance/rejectance rule. 

The function {\tt XLALMCMCJump} computes a proposal parameter set from the current parameter set, using the covariance matrix and a random draw from a multidimensional Student distribution:
The use of a Student distribution (with n=2) ensures that outliers are more weighted than in a Normal distribution.

The function {\tt XLALMCMCCheckAnnealing} is used for annealing if the flag {\tt useAnnealing} is set to 1. In this function, an 'annealing temparature' {\tt annealingTemp} decreases from an initial temperature of {\tt annealingTempBegin} to unity within {\tt annealingSteps} stepsin an exponential way. 

The function {\tt XLALMCMCCheckBurnin} is experimental code to check if the burnin is reached, if {\tt flagBurnin} is set to 1. The last {\tt burninStep} values of the chain are fitted to a line every {\tt burninNumber} steps in the chain, and if the ratio of slope to mean is below a value {\tt burninThreshold} the burnin is reached.  

The function {\tt XLALMCMCCheckUpdate} is used to update the covariance matrix. The updating is done by using the actual parameters and a set of mean values ({\tt mean}) and difference values ({\tt xdiff} and {\tt ndiff} ). 

The function {\tt XLALMCMCCheckAcceptRatio} is experimental code to check the current acceptance ratio. If this ratio is too small. If, after {\tt acceptRatioNorm} steps the ratio of accepted ({\tt acceptRatioCounter}) to rejected jumps is below {\tt acceptThreshold}, then the annealing temperature is halfed, and the counters are zeroed.

The functions {\tt XLALMCMCAddParam, XLALMCMCGetParam, XLALMCMCGetParameter, XLALMCMCSetParameter, XLALMCMCCopyPara, XLALMCMCFreePara} and {\tt XLALMXMXDestroyPara} are functions to handle the {\tt MCMCParameter} structure. They can be used to add new parameter to this structure (in the initializing step ({\tt XLALMCMCAddParam}) and to set/get parameters. The function {\tt XLALMCMCCopyPara} copies the parameter structure (actually only the current values. The pointer to the core-structure remains the same).

The functions {\tt XLALMultiStudentDeviates} and {\tt XLALMultiNormalDeviates} draw multivariate random values, either Student-distributed or Norla distributed, given a covariance matrix. A function {\tt XLALCheckPositiveDefinite} can be used to check if the function is positive definite, which is required beforehand.


To set a different chain for the parameter estimation, just set another random seed. 

\subsubsection*{Algorithm}

The algorithms used in these functions are explained in detail in [Ref Needed].

\subsubsection*{Uses}

This section briefly explains how to set-up the \texttt{MCMCInput} and the \texttt{MCMCParameter} structure for doing a MCMC sampling. 

The first step is to initialize the \texttt{MCMCInput} structure with the output from data reading and conditioning, as well with user inputs to define settings for the MCMC. Also, the user is required to set the pointer to three functions which will initialize the parameter structure ({\tt funcInit}), that calculates the logarithm of the likelihood for a given set of parameters ({\tt funcLikelihood}) and that calculates the logarithm of the prior given a set of parameters ({\tt funcPrior}). Below is an example:

\begin{verbatim}
        LALMCMCInput inputMCMC;
        inputMCMC.fcFilterInput  = fcFilterInput;
        inputMCMC.fcFilterParams = fcFilterParams;
        inputMCMC.fcTmpltParams  = fcTmpltParams;
        inputMCMC.fcDataParams   = fcDataParams;
        inputMCMC.randParams     = randParams;
        inputMCMC.approximant    = approximant;
        inputMCMC.inspiralTable  = inputCurrent;

        inputMCMC.tmpltPtr = 
          (InspiralTemplate*)LALCalloc(sizeof(InspiralTemplate),1);

        inputMCMC.counter=0;

        inputMCMC.useAnnealing = 0;        
        inputMCMC.numberAnneal=iterAnneal;
        inputMCMC.annealingTempBegin=4.0;
        inputMCMC.annealingTemp=1.0;
        inputMCMC.annealingSteps=0;

        inputMCMC.useScaling = 1;        
        inputMCMC.scalePeak=50.0;
        inputMCMC.scaleNormal=1.0;
        inputMCMC.scaleQ=0.84;
        inputMCMC.scalePA=1e-3;

        inputMCMC.flagAcceptRatio=1;
        inputMCMC.acceptRatioCounter=0;

        inputMCMC.useUpdate=1;
        inputMCMC.updateNumber=iterCovupdate;
        inputMCMC.updateOffset=100;
        inputMCMC.updateCounter=0;
        inputMCMC.mean=NULL;
        inputMCMC.xdiff=NULL;
        inputMCMC.ndiff=NULL;

        inputMCMC.flagBurnin = 1;        
        inputMCMC.burninNumber = 100;
        inputMCMC.burninStep = 100;
        inputMCMC.burninTime = 0;

        inputMCMC.numberDraw=iterDraw;

        inputMCMC.verbose=1;

        inputMCMC.funcInit     = MCMCInit0;
        inputMCMC.funcTemplate = MCMCTemplate0;
        inputMCMC.funcPrior    = MCMCPrior0;

\end{verbatim}



The next step is to initialize a \texttt{MCMCParameter} structure as a NULL pointer and pass it along the \texttt{MCMCInput} structure to \texttt{XLALMCMCMetro}. Thats it!

However, you need to point to a function that initializes the parameter structure. This \texttt{funcInit} function populates the parameter structure with as many parameters as the user whishes, defining inital and boundary values. An example of a init function is given here, which also can be found in the code {\tt XLALMCMCUser.c}:

\subsubsection*{Notes}


</lalLaTeX>  */
#endif

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include "LALInspiralMCMC.h"
#include <lal/Date.h>
#include <lal/Random.h>
#include <lal/AVFactories.h>
#include <lal/TimeDelay.h>
#include <lal/DetResponse.h>


#define rint(x) floor((x)+0.5)

NRCSID (LALINSPIRALMCMCC, "$Id: LALInspiralPhase.c,v 1.9 2003/04/14 00:27:22 sathya Exp $"); 

/* *****************************
printMatrix
  ***************************** */
void 
printMatrix( gsl_matrix *covMat, int dim)
{
  /* print matrix content */
  int i, j;

  for (i = 0; i < dim; i++) {
    printf("MCMCINFO: ");
    for (j = 0; j < dim; j++) {
      printf("  %10e", gsl_matrix_get ( covMat, i,j));
    }
    printf("\n");
  }
}

/* *****************************
printState
  ***************************** */
void 
printState( StateMCMC state )
{
  switch (state)
  {
  case unknownState:
    printf("unknown state");
    break;
  case doScaling:
    printf("scaling");
    break;
  case doAnnealing:
    printf("annealing");
    break;
  case doUpdating:
    printf("maxtrix updating");
    break;
  case doDrawing:
    printf("draw values");
    break;
  }
}



/* *****************************
XLALInspiralMCMCMetro
  ***************************** */
/*  <lalVerbatim file="XLALMCMCMetroCP"> */
void
XLALMCMCMetro ( 
  LALMCMCParameter  **paraPtr,
  LALMCMCInput       *inputMCMC
) { /* </lalVerbatim>  */

  static LALStatus status;

  LALMCMCParameter *parameter   = NULL;
  LALMCMCParam     *paraHead    = NULL;
  gsl_matrix       *startCovMat = NULL;
  gsl_matrix       *covMat      = NULL;
  REAL4Vector      *vector      = NULL;
  StateMCMC state;
  REAL4     logPrior, logLikelihood, currentLogPosterior;
  REAL8     range=0.0;
  int       foundBurnin, loopMCMC, move, dim, i,j;

  /* initialize the parameter structure */
  parameter = *paraPtr;  
  inputMCMC->funcInit( parameter, inputMCMC->inspiralTable ); 

  dim            = parameter->dimension;
  inputMCMC->dim = dim;  

  /* initialize the internal vectors */
  for (paraHead=parameter->param; paraHead; 
       paraHead=paraHead->next) 
  {
    paraHead->core->chain=NULL;
    LALSCreateVector( &status, &paraHead->core->chain, inputMCMC->numberDraw+1);
  }                
 
  /* initialize vectors in the MCMC structure */
  LALSCreateVector( &status, &vector, dim);
  inputMCMC->mean = (REAL8*)LALCalloc( sizeof(REAL8), dim );
  inputMCMC->xdiff= (REAL8*)LALCalloc( sizeof(REAL8), dim );
  inputMCMC->ndiff= (REAL8*)LALCalloc( sizeof(REAL8), dim );

  /* prepare the gsl matrices */
  startCovMat = gsl_matrix_alloc( dim, dim); 
  covMat      = gsl_matrix_alloc( dim, dim); 

  /* set everything to zeros in the matrix */
  for (i=0;i<dim;i++) 
    for (j=0;j<dim;j++) 
      gsl_matrix_set( startCovMat, i, j, 0.0 );

  /* populate the startCovMat with initial values */
  printf("populating the covariance matrix:\n");
  for (paraHead=parameter->param,i=0; paraHead; paraHead=paraHead->next,i++)
  {    
    range=paraHead->core->maxVal - paraHead->core->minVal;
    gsl_matrix_set( startCovMat, i, i, 0.001*range*range );
    printf("element %d: %f\n", i, 0.001*range*range);
  }
  gsl_matrix_memcpy( covMat, startCovMat ); 

  /* initialize the state */
  inputMCMC->counter=inputMCMC->counterState=inputMCMC->counterAccept=inputMCMC->counterAcceptDraw=0;
  inputMCMC->burninCounter=0;
  loopMCMC=1;
  if (inputMCMC->useScaling)
  {
    state = doScaling;
    inputMCMC->scaling = inputMCMC->scalePeak;
  }
  else if (inputMCMC->useAnnealing)
  {
    state = doAnnealing;
    inputMCMC->scaling = inputMCMC->scaleNormal;
  }
  else
  {
    state = unknownState; /* TODO: error*/
    printf("ERROR: unknown state\n");
    exit(1);
  }
      

  /* get the first values */
  inputMCMC->funcPrior( inputMCMC, parameter );
  logLikelihood = inputMCMC->funcLikelihood( inputMCMC, parameter);
  logPrior = parameter->logPrior;
  currentLogPosterior=logLikelihood+logPrior;

  do
  {
    /* increase the counters */
    inputMCMC->counter++;      /*k*/
    inputMCMC->counterState++; /*c*/

    /* print generel information */
    printf("####################################################\n");
    printf("\nMCMCSTATUS k: %d c: %d  state: ", 
           inputMCMC->counter, inputMCMC->counterState);
    printState( state);
    printf("\n");
    
   
    /* do the sampling from the underlying distribution */
    move=XLALMCMCSample(inputMCMC, &parameter, &currentLogPosterior, covMat);

    /* count the accepted moves */
    if (move) 
    {
      inputMCMC->counterAccept++;
    }
    

    /* ------------------------------
       checking the burnin
       ------------------------------ */
    if ( inputMCMC->flagBurnin  && state == doScaling)
    {

      /* store the actual parameters in the chain structure */
      for (paraHead=parameter->param; paraHead; paraHead=paraHead->next)
      {
        paraHead->core->chain->data[inputMCMC->burninCounter]= paraHead->value;
      }

      /* increase internal burnin counter */
      inputMCMC->burninCounter++;

      if ( inputMCMC->counterState>=inputMCMC->burninNumber && 
           inputMCMC->burninCounter>=inputMCMC->burninStep )
      {
        /* reset internal burnin counter */
        inputMCMC->burninCounter=0;

        /* check if burnin found */
        foundBurnin=XLALMCMCCheckBurnin( inputMCMC, parameter );  

        if (foundBurnin || inputMCMC->counterState >= inputMCMC->burninMaxNumber ) 
        {
          if (inputMCMC->verbose)
            printf("MCMCINFO: Burnin period found at step %d." 
                   "Moving to updating state.\n", inputMCMC->counter );
          
          /* reset the counter */
          inputMCMC->burninTime = inputMCMC->counter;
          
          /* adjust scaling factor and the cov-matrix*/
          inputMCMC->scaling = inputMCMC->scaleNormal;
          gsl_matrix_scale( covMat, 0.20);
          
          /* deactivate this algorithm and set the new state */
          inputMCMC->flagBurnin=0;          
          state=doUpdating;
          inputMCMC->counterState=0;          
        } 
        
      }
    }

    /* print content of the covariance matrix to screen
       if verbose flag is set */
    if (inputMCMC->verbose)
    {
      printMatrix( covMat, dim );
    }

    /* ------------------------------
       updating the covariance matrix 
       ------------------------------ */
    if ( state == doUpdating )
    {
      /* check if enough iterations used here */
      if ( inputMCMC->counterState < inputMCMC->updateNumber )
      {
        /* update the covariance matrix */
        XLALMCMCCheckUpdate( inputMCMC, parameter, covMat );
      }
      else
      {
        /* move to the next state */
        if (inputMCMC->verbose)
        {
          printf("MCMCINFO: Matrix updating complete at step %d. "
                 "Moving to drawing state. "                        \
                 "Content of the matrix: \n", inputMCMC->counter);
          printMatrix(covMat, dim);
        }

        /* set drawing state after matrix updating is finished */
        state = doDrawing;       
        inputMCMC->counterState=0;
      }

    }

    /* ------------------------------
       drawing the values from the chain, and store them
       ------------------------------ */
    if ( state == doDrawing )
    {
      /* store the drawn values in the chain */
      for (paraHead=parameter->param; paraHead; paraHead=paraHead->next)
      {
        paraHead->core->chain->data[inputMCMC->counterState]= paraHead->value;
      }

      /* check if end of loop is reached, i.e. if enough values has been drawn */
      if ( inputMCMC->counterState == inputMCMC->numberDraw )
      {
        loopMCMC=0; /* ending condition for this loop */
      }
    }


    /*if (inputMCMC->counterState>1000) loopMCMC=0;*/
    
  } while( loopMCMC );


  LALFree(inputMCMC->mean);
  LALFree(inputMCMC->xdiff);
  LALFree(inputMCMC->ndiff);
  gsl_matrix_free(covMat);
  gsl_matrix_free( startCovMat);

  /* print some summary informations */
  if ( inputMCMC->verbose>0 )
  {
    printf("MCMCINFO: Total number of iterations: %d , steps accepted: %d  ( %.1f %% )\n", 
           inputMCMC->counter, inputMCMC->counterAccept, 
           100.0*(float)inputMCMC->counterAccept/(float)inputMCMC->counter );
    printf("MCMCINFO: Number of draws used: %d, steps accepted: %d ( %.1f %%)\n",
           inputMCMC->numberDraw, inputMCMC->counterAccept, 
           100.0*(float)inputMCMC->counterAccept/(float)inputMCMC->numberDraw );    
  }

  /* set the correct pointer */
  paraPtr= &parameter;
}

/* *****************************
  XLALMCMCCheckAnnealing
  ***************************** */
/*  <lalVerbatim file="XLALMCMCCheckAnnealingCP"> */
void 
XLALMCMCCheckAnnealing(
  gsl_matrix *covMat,
  gsl_matrix *inputMat,
  LALMCMCInput *inputMCMC) 
{ /* </lalVerbatim>  */

  int k, number;
  k=inputMCMC->counterState;
  number = inputMCMC->annealingSteps;

  /* adjust the 'temperature' of the MCMC every annealingSteps steps */
  if ( (k%number)==0 ) 
  {
    /* calculate the 'heat' */
    inputMCMC->annealingTemp = exp( log(inputMCMC->annealingTempBegin) * (number-k)/(number-1.0));
    
    /* copy and scale matrix */
    gsl_matrix_memcpy( covMat, inputMat );
    gsl_matrix_scale( covMat, inputMCMC->annealingTemp );
    
  }  
}


/* *****************************
  XLALMCMCCheckBurnin
  ***************************** */
/*  <lalVerbatim file="XLALMCMCCheckBurninCP"> */
INT4
XLALMCMCCheckBurnin(
  LALMCMCInput *inputMCMC,
  LALMCMCParameter *parameter)
{ /* </lalVerbatim>  */

  LALMCMCParam *paraHead=NULL;

  UINT4 i,n;
  REAL4 x,y, sxy, sxx, sx, sy;
  REAL4 mean, slope;
  INT4 flag;


  flag=1;
  n=inputMCMC->burninStep;


  printf("\nXLALMCMCCheckBurnin: -------------------------------------------\n");
  /* get the data for each parameter and perform a linear regression */
  for (paraHead=parameter->param; paraHead; paraHead=paraHead->next)
  { 

    sxy=0;
    sxx=0;
    sx=0;
    sy=0;

    for (i=0; i<n;i++) {
      x=(float)i;
      y=paraHead->core->chain->data[i];
      sxy+=x*y;
      sxx+=x*x;
      sx+=x;
      sy+=y;
    }

    /* calculate the slope and the mean */
    mean=sy/(float)n;
    slope=((float)n*sxy-sx*sy)/((float)n*sxx-sx*sx);
    
    /* is fluctuation is above a trheshold, set flag to zero: no burnin */
    if (slope/mean>inputMCMC->burninThreshold) 
    {
      flag=flag*0;
    } 

    printf("  parameter %10s: mean=%10e  slope=%10e \n", paraHead->core->name, mean, slope);

  }
  printf("-------------------------------------------\n\n");

  return flag;
}

/* *****************************
  XLALMCMCCheckUpdate
  ***************************** */
/*  <lalVerbatim file="XLALMCMCCheckUpdateCP"> */
void
XLALMCMCCheckUpdate(
   LALMCMCInput *inputMCMC,
   LALMCMCParameter *parameter,
   gsl_matrix *covMat)
{ /* </lalVerbatim>  */

  LALMCMCParam* paraHead=NULL;
  REAL8* mean = NULL;
  REAL8* oldDiff= NULL;
  REAL8* newDiff = NULL;
  UINT4 i,j,c;
  REAL8 element, term1, term2, t;
  REAL8 minWrap, maxWrap, deltaWrap;

  /* set some frequent used parameters */
  mean=inputMCMC->mean;
  oldDiff=inputMCMC->xdiff;
  newDiff=inputMCMC->ndiff;

  /* increase internal counter */
  inputMCMC->updateCounter++;
  
  /* calculate oldDiff */ 
  for (paraHead=parameter->param,c=0; paraHead; paraHead=paraHead->next,c++)
  {

    if ( inputMCMC->updateCounter == 1)
    {
      mean[c] = paraHead->value;
      oldDiff[c]=0.0;
      newDiff[c]=0.0;
    }
    else      
    {
      /* calculate oldDiff */
      oldDiff[c] = paraHead->value - mean[c];
      
      /* check wrapping on the oldDiff value*/        
      minWrap=paraHead->core->minVal;
      maxWrap=paraHead->core->maxVal;
      deltaWrap=maxWrap-minWrap;
      if (paraHead->core->wrapping)
      {
        if (oldDiff[c] > maxWrap) oldDiff[c]-= deltaWrap;
        if (oldDiff[c] < minWrap) oldDiff[c]+= deltaWrap;
      }
      
      /* update mean */
      mean[c] += oldDiff[c]/(inputMCMC->updateCounter + inputMCMC->updateOffset);
      
      /* check wrapping on the mean value */
      if (paraHead->core->wrapping)
      {
        if (mean[c] > maxWrap) mean[c]-= deltaWrap;
        if (mean[c] < minWrap) mean[c]+= deltaWrap;
      }		    
      
      /* calculate newDiff (usage of NEW mean) */
      newDiff[c] = paraHead->value - mean[c];
      if (paraHead->core->wrapping)
      {
        if (newDiff[c] > maxWrap) newDiff[c]-= deltaWrap;
        if (newDiff[c] < minWrap) newDiff[c]+= deltaWrap;
      }
    }
  }
     
  /* update covariance: */ 
  t=(float)inputMCMC->updateCounter-1.0;
  for (i=0; i<inputMCMC->dim; i++)
    for (j=0; j<inputMCMC->dim; j++)
      if (t>1) 
      {              
        element=(t-1)*gsl_matrix_get( covMat, i, j);
        term1 = t*oldDiff[i]*oldDiff[j]/((t+1.0)*(t+1.0));
        term2 = newDiff[i]*newDiff[j];            
        gsl_matrix_set( covMat, i,j, (element+term1+term2)/t );	  
      }      
     
}


/* ******************************************
  XLALMCMCSample
  ******************************************* */
/*  <lalVerbatim file="XLALMCMCSampleCP"> */
UINT4 
XLALMCMCSample(
  LALMCMCInput *inputMCMC,
  LALMCMCParameter **paraPtr,
  REAL4  *oldLogPosterior, 
  gsl_matrix *covMat
  )
{ /* </lalVerbatim>  */

  static LALStatus status;

  LALMCMCParameter *parameter=NULL; /* parameter value    */
  LALMCMCParameter *proposal=NULL;  /* proposal  value    */
  LALMCMCParameter *help=NULL;      /* help parameter set */
  LALMCMCParam* paraHead = NULL;

  REAL4 alpha, my_random, s; 
  REAL4 logPrior, logLikelihood, logPosterior;
  UINT4 move, accept, c;
  INT4 testPrior;

  /* set the parameter */
  parameter=*paraPtr;
  move=0;

  /* allocate spec for proposal set */
  proposal=(LALMCMCParameter*)LALMalloc( sizeof(LALMCMCParameter) ); 

  do {
    accept = 1;
    XLALMCMCCopyPara( &proposal, parameter);
    XLALMCMCJump( inputMCMC, proposal, covMat); 
    
    for (paraHead = proposal->param,c=0; paraHead; paraHead=paraHead->next,c++ )
    {
      /* check if parameter lies in valid range */
      if (paraHead->value < paraHead->core->minVal || 
          paraHead->value > paraHead->core->maxVal)
      {
        accept=0;
        /*
        if ( inputMCMC->verbose )
          printf("MCMCSAMPLE Parameter %10s outside range. Value: %8.3f  Range: %8.3f - %8.3f \n", 
                 paraHead->core->name, paraHead->value,
                 paraHead->core->minVal,  paraHead->core->maxVal );
        */
      }
    }
  } while ( !accept );

  /* calculate the log Prior */
  testPrior = inputMCMC->funcPrior( inputMCMC, proposal );
  logPrior = proposal->logPrior;
  
  /*-- determine likelihood if prior gives something reasonable:   --*/
  if ( logPrior>-HUGE_VAL ) 
  {
    
    /* calculate the new posterior value for the proposal parameter set */
    logLikelihood = inputMCMC->funcLikelihood( inputMCMC, proposal );
    logPosterior = logLikelihood + logPrior;
    
    /* calculate the alpha-value and draw a random number */
    s=inputMCMC->scaling;
    alpha=exp( s*(logPosterior - *oldLogPosterior) );
    LALUniformDeviate( &status, &my_random, inputMCMC->randParams );
    
    /* now check accept/reject criterion */
    if ( my_random<=alpha ) 
    {        
      /* accept the proposal set */
      move = 1;
      
      /* just swap the two pointers */
      help=proposal;
      proposal=parameter;
      parameter=help;

      /* return the new log posterior value */
      *oldLogPosterior = logPosterior;
    } 

    /* output */ 
    if ( inputMCMC->verbose )
    {
      if (move==1)
        printf("MCMCSAMPLE: ++JumpIsAccepted ");
      else
        printf("MCMCSAMPLE: --JumpNotAccepted ");
      printf("current logPost: %6.3f  proposal logPost: %6.3f alpha: %6.3f  "
             "u: %6.3f\n", 
             *oldLogPosterior, logPosterior, alpha, my_random);
    }
  }
 
  /* test printout */  
  if ( inputMCMC->verbose )
  {
    printf("MCMCPARAMETER: ");
    printf("| SNR: %f  ", *oldLogPosterior);  
    for (paraHead=parameter->param; paraHead; paraHead=paraHead->next) {
      printf(" | %s: %9.5f", paraHead->core->name, paraHead->value);
    }
    fprintf(stdout, "\n");
  }
  

  /* recopy the correct parameter structure */
  *paraPtr=parameter;

  /* free proposal parameter set */
  XLALMCMCFreePara( proposal );

  return move;
}

void XLALMCMCGetCartesianPos(REAL8 vec[3],LALMCMCParameter *parameter)
{
REAL8 longitude = XLALMCMCGetParameter(parameter,"long");
REAL8 latitude = XLALMCMCGetParameter(parameter,"lat");
/*REAL8 distance = XLALMCMCGetParameter(parameter,"distMpc");*/
vec[0]=cos(longitude)*cos(latitude);
vec[1]=sin(longitude)*cos(latitude);
vec[1]=sin(latitude);
return;
}

void CartesianToSkyPos(REAL8 pos[3],LALMCMCParameter *parameter)
{
REAL8 longi,lat,dist;
dist=sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
/*XLALMCMCSetParameter(parameter,"distMpc",dist);*/
longi=atan2(pos[1]/dist,pos[0]/dist);
	if(longi<0.0) longi=LAL_TWOPI+longi;
lat=asin(pos[2]/dist);
XLALMCMCSetParameter(parameter,"lat",lat);
XLALMCMCSetParameter(parameter,"long",longi);
return;
}

void crossProduct(REAL8 out[3],REAL8 x[3],REAL8 y[3])
{
out[0]=x[1]*y[2] - x[2]*y[1];
out[1]=y[0]*x[2] - x[0]*y[2];
out[2]=x[0]*y[1] - x[1]*y[0];
return;
}

void normalise(REAL8 vec[3]);
void normalise(REAL8 vec[3]){
REAL8 my_abs=0.0;
my_abs=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
vec[0]/=my_abs;
vec[1]/=my_abs;
vec[2]/=my_abs;
return;
}

INT4 XLALMCMCDifferentialEvolution(
	LALMCMCInput *inputMCMC,
	LALMCMCParameter *parameter)
{
	static LALStatus status;
	LALMCMCParameter **Live=inputMCMC->Live;
	int i=0,j=0,dim=0,same=1;
	REAL4 randnum;
	int Nlive = (int)inputMCMC->Nlive;
	LALMCMCParam *paraHead=NULL;
	LALMCMCParam *paraA=NULL;
	LALMCMCParam *paraB=NULL;

	dim = parameter->dimension;
	if(inputMCMC->randParams==NULL) LALCreateRandomParams(&status,&(inputMCMC->randParams),0);
	/* Select two other samples A and B*/
	LALUniformDeviate(&status,&randnum,inputMCMC->randParams);
	i=(int)(Nlive*randnum);
	/* Draw two different samples from the basket. Will loop back here if the original sample is chosen*/
	drawtwo:
	do {LALUniformDeviate(&status,&randnum,inputMCMC->randParams); j=(int)(Nlive*randnum);} while(j==i);
	paraHead=parameter->param;
	paraA=Live[i]->param; paraB=Live[j]->param;
		/* Add the vector B-A */
	same=1;
	while(paraHead)
	{
		if(paraHead->value!=paraA->value && paraHead->value!=paraB->value) same=0;
		paraHead->value+=paraB->value-paraA->value;
		paraB=paraB->next; paraA=paraA->next;
		paraHead=paraHead->next;
	}
	if(same==1) goto drawtwo;
	/* Bring the sample back into bounds */
	XLALMCMCCyclicReflectiveBound(parameter);
	return(0);
}

INT4 XLALMCMCReflectDetPlane(
	LALMCMCInput *inputMCMC,
	LALMCMCParameter *parameter
	)
{ /* Function to reflect a point on the sky about the plane of 3 detectors */
  /* Returns -1 if not possible */
static LALStatus status;
UINT4 i,j;
int DetCollision=0;
REAL4 randnum;
REAL8 longi,lat,newlong,newlat;
REAL8 dist;
REAL8 pos[3];
REAL8 normal[3];
REAL8 w1[3]; /* work vectors */
REAL8 w2[3];
INT4 IFO1,IFO2,IFO3;
REAL8 detvec[3];

if(inputMCMC->numberDataStreams<3) return(-1) ; /* not enough IFOs to construct a plane */
for(i=0;i<inputMCMC->numberDataStreams;i++)
	for(j=i;j<inputMCMC->numberDataStreams;j++)
		if((j!=i) && inputMCMC->detector[i]==inputMCMC->detector[j]) DetCollision+=1;

if(inputMCMC->numberDataStreams-DetCollision<3) return(-1); /* Not enough independent IFOs */

/* Select IFOs to use */
if(inputMCMC->randParams==NULL) LALCreateRandomParams(&status,&(inputMCMC->randParams),0);
LALUniformDeviate(&status,&randnum,inputMCMC->randParams);
IFO1 = (INT4)floor(inputMCMC->numberDataStreams*randnum);	
LALUniformDeviate(&status,&randnum,inputMCMC->randParams);
IFO2 = (INT4)floor((inputMCMC->numberDataStreams-1)*randnum);
while(IFO1==IFO2 || inputMCMC->detector[IFO1]==inputMCMC->detector[IFO2]) IFO2=(IFO2+1) % inputMCMC->numberDataStreams;
LALUniformDeviate(&status,&randnum,inputMCMC->randParams);
IFO3 = (INT4)floor((inputMCMC->numberDataStreams-2)*randnum);
while(IFO3==IFO1
	|| IFO3==IFO2
	|| inputMCMC->detector[IFO3]==inputMCMC->detector[IFO1]
	|| inputMCMC->detector[IFO3]==inputMCMC->detector[IFO2])
	IFO3=(IFO3+1) % inputMCMC->numberDataStreams;
/*fprintf(stderr,"Using %s, %s and %s for plane\n",inputMCMC->ifoID[IFO1],inputMCMC->ifoID[IFO2],inputMCMC->ifoID[IFO3]);*/

longi = XLALMCMCGetParameter(parameter,"long");
lat = XLALMCMCGetParameter(parameter,"lat");

double deltalong=0;

/* Convert to earth coordinates */
	SkyPosition geodetic,equatorial;
	equatorial.longitude=longi;
	equatorial.latitude=lat;
	equatorial.system=COORDINATESYSTEM_EQUATORIAL;
	geodetic.system=COORDINATESYSTEM_GEOGRAPHIC;
	LALEquatorialToGeographic(&status,&geodetic,&equatorial,&(inputMCMC->epoch));
	deltalong=geodetic.longitude-equatorial.longitude;

XLALMCMCSetParameter(parameter,"long",deltalong+XLALMCMCGetParameter(parameter,"long"));
XLALMCMCGetCartesianPos(pos,parameter); /* Get sky position in cartesian coords */


/* calculate the unit normal vector of the detector plane */
for(i=0;i<3;i++){ /* Two vectors in the plane */
	w1[i]=inputMCMC->detector[IFO2]->location[i] - inputMCMC->detector[IFO1]->location[i];
	w2[i]=inputMCMC->detector[IFO3]->location[i] - inputMCMC->detector[IFO1]->location[i];
	detvec[i]=inputMCMC->detector[IFO1]->location[i];
	}
crossProduct(normal,w1,w2);
normalise(normal);
normalise(detvec);

/* Calculate the distance between the point and the plane n.(point-IFO1) */
for(dist=0.0,i=0;i<3;i++) dist+=pow(normal[i]*(pos[i]-detvec[i]),2.0);
dist=sqrt(dist);
/* Reflect the point pos across the plane */
for(i=0;i<3;i++) pos[i]=pos[i]-2.0*dist*normal[i];


CartesianToSkyPos(pos,parameter);
XLALMCMCSetParameter(parameter,"long",XLALMCMCGetParameter(parameter,"long")-deltalong);

	/* Compute change in tgeocentre for this change in sky location */
	newlong=XLALMCMCGetParameter(parameter,"long");
	newlat=XLALMCMCGetParameter(parameter,"lat");
	REAL8 dtold,dtnew,deltat;
	DetTimeAndASource DTAAS; /* This holds the source and the detector */
	LALSource source; /* The position and polarisation of the binary */
	source.equatorialCoords.longitude = longi;
	source.equatorialCoords.latitude = lat;
	source.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
	LALPlaceAndGPS det_gps; /* This will hold the detector site and epoch of observation */
	det_gps.p_gps=&(inputMCMC->epoch);
	DTAAS.p_source = &(source.equatorialCoords);
	DTAAS.p_det_and_time=&det_gps;
	DTAAS.p_det_and_time->p_detector = inputMCMC->detector[0]; /* Select detector */
	LALTimeDelayFromEarthCenter(&status,&dtold,&DTAAS); /* Compute time delay */
	source.equatorialCoords.longitude=newlong;
	source.equatorialCoords.latitude=newlat;
	LALTimeDelayFromEarthCenter(&status,&dtnew,&DTAAS); /* Compute time delay */
	deltat=dtold-dtnew; /* deltat is change in arrival time at geocentre */
	deltat+=XLALMCMCGetParameter(parameter,"time");
	XLALMCMCSetParameter(parameter,"time",deltat);

XLALMCMCCyclicReflectiveBound(parameter);

return(0);
}

void XLALMCMCRotateSky(
	LALMCMCInput *inputMCMC,
	LALMCMCParameter *parameter
	)
{ /* Function to rotate the current sample around the vector between two random detectors */
	static LALStatus status;
	INT4 IFO1,IFO2;
	REAL4 randnum;
	REAL8 vec[3];
	REAL8 cur[3];
	REAL8 longi,lat;
	REAL8 vec_abs=0.0,theta,c,s;
	INT4 i,j;

	if(inputMCMC->numberDataStreams<2) return;
	if(inputMCMC->numberDataStreams==2 && inputMCMC->detector[0]==inputMCMC->detector[1]) return;
	
	longi = XLALMCMCGetParameter(parameter,"long");
	lat = XLALMCMCGetParameter(parameter,"lat");
	
	/* Convert the RA/dec to geodetic coordinates, as the detectors use these */
	SkyPosition geodetic,equatorial;
	equatorial.longitude=longi;
	equatorial.latitude=lat;
	equatorial.system=COORDINATESYSTEM_EQUATORIAL;
	geodetic.system=COORDINATESYSTEM_GEOGRAPHIC;
	LALEquatorialToGeographic(&status,&geodetic,&equatorial,&(inputMCMC->epoch));
	longi=geodetic.longitude;
	lat=geodetic.latitude;
	cur[0]=cos(lat)*cos(longi);
	cur[1]=cos(lat)*sin(longi);
	cur[2]=sin(lat);
	
	if(inputMCMC->randParams==NULL) LALCreateRandomParams(&status,&(inputMCMC->randParams),0);
	LALUniformDeviate(&status,&randnum,inputMCMC->randParams);
	IFO1 = (INT4)floor(inputMCMC->numberDataStreams*randnum);
	do{ /* Pick random interferometer other than the first one */
		LALUniformDeviate(&status,&randnum,inputMCMC->randParams);
		IFO2 = (INT4)floor(inputMCMC->numberDataStreams*randnum);
	}while(IFO2==IFO1 || inputMCMC->detector[IFO1]==inputMCMC->detector[IFO2]);
	
/*	fprintf(stderr,"Rotating around %s-%s vector\n",inputMCMC->ifoID[IFO1],inputMCMC->ifoID[IFO2]);*/
	/* Calc normalised direction vector */
	for(i=0;i<3;i++) vec[i]=inputMCMC->detector[IFO2]->location[i]-inputMCMC->detector[IFO1]->location[i];
	for(i=0;i<3;i++) vec_abs+=vec[i]*vec[i];
	vec_abs=sqrt(vec_abs);
	for(i=0;i<3;i++) vec[i]/=vec_abs;
	
	/* Chose random rotation angle */
	LALUniformDeviate(&status,&randnum,inputMCMC->randParams);
	theta=LAL_TWOPI*randnum;
	c=cos(-theta); s=sin(-theta);
	/* Set up rotation matrix */
	double R[3][3] = {{c+vec[0]*vec[0]*(1.0-c), 
                     vec[0]*vec[1]*(1.0-c)-vec[2]*s,
                     vec[0]*vec[2]*(1.0-c)+vec[1]*s},
                    {vec[1]*vec[0]*(1.0-c)+vec[2]*s,
                     c+vec[1]*vec[1]*(1.0-c),
                     vec[1]*vec[2]*(1.0-c)-vec[0]*s},
                    {vec[2]*vec[0]*(1.0-c)-vec[1]*s,
                     vec[2]*vec[1]*(1.0-c)+vec[0]*s,
                     c+vec[2]*vec[2]*(1.0-c)}};
	REAL8 new[3]={0.0,0.0,0.0};
	for (i=0; i<3; ++i)
		for (j=0; j<3; ++j)
			new[i] += R[i][j]*cur[j];
	double newlong = atan2(new[1],new[0]);
	if(newlong<0.0) newlong=LAL_TWOPI+newlong;
	
	geodetic.longitude=newlong;
	geodetic.latitude=asin(new[2]);
	/* Convert back into equatorial (sky) coordinates */
	LALGeographicToEquatorial(&status,&equatorial,&geodetic,&(inputMCMC->epoch));
	newlong=equatorial.longitude;
	double newlat=equatorial.latitude;
	
	/* Compute change in tgeocentre for this change in sky location */
	REAL8 dtold,dtnew,deltat;
	DetTimeAndASource DTAAS; /* This holds the source and the detector */
	LALSource source; /* The position and polarisation of the binary */
	source.equatorialCoords.longitude = longi;
	source.equatorialCoords.latitude = lat;
	source.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
	LALPlaceAndGPS det_gps; /* This will hold the detector site and epoch of observation */
	det_gps.p_gps=&(inputMCMC->epoch);
	DTAAS.p_source = &(source.equatorialCoords);
	DTAAS.p_det_and_time=&det_gps;
	DTAAS.p_det_and_time->p_detector = inputMCMC->detector[0]; /* Select detector */
	LALTimeDelayFromEarthCenter(&status,&dtold,&DTAAS); /* Compute time delay */
	source.equatorialCoords.longitude=newlong;
	source.equatorialCoords.latitude=newlat;
	LALTimeDelayFromEarthCenter(&status,&dtnew,&DTAAS); /* Compute time delay */
	deltat=dtold-dtnew; /* deltat is change in arrival time at geocentre */
	deltat+=XLALMCMCGetParameter(parameter,"time");
	XLALMCMCSetParameter(parameter,"time",deltat);	
	XLALMCMCSetParameter(parameter,"lat",newlat);
	XLALMCMCSetParameter(parameter,"long",newlong);
	/*fprintf(stderr,"Skyrotate: new pos = %lf %lf %lf => %lf %lf\n",new[0],new[1],new[2],newlong,asin(new[2]));*/
	XLALMCMCCyclicReflectiveBound(parameter);

	return;
}

int XLALMCMC1PNMasseta(LALMCMCInput *inputMCMC, LALMCMCParameter *parameter)
{
	REAL8 eta1,eta2,mc1,mc2;
	REAL4 randnum;
	int logflag=0;
	static LALStatus status;
	eta1=XLALMCMCGetParameter(parameter,"eta");
	if(XLALMCMCCheckParameter(parameter,"logM")) {mc1=exp(XLALMCMCGetParameter(parameter,"logM")); logflag=1;}
	else mc1=XLALMCMCGetParameter(parameter,"mchirp");
	LALUniformDeviate(&status,&randnum,inputMCMC->randParams);
	eta2=0.25*(REAL8)randnum;
	mc2 = pow(eta2/eta1,3./5.)*mc1;
	XLALMCMCSetParameter(parameter,"eta",eta2);
	if(logflag) XLALMCMCSetParameter(parameter,"logM",log(mc2));
	else XLALMCMCSetParameter(parameter,"mchirp",mc2);
	
	return(0);
}


/* ******************************************
  XLALMCMCJump
  ******************************************* */
/*  <lalVerbatim file="XLALMCMCJumpCP"> */
void 
XLALMCMCJump(
  LALMCMCInput     *inputMCMC,
  LALMCMCParameter *parameter, 
  gsl_matrix       *covMat
  ) 
{ /* </lalVerbatim>  */
  static LALStatus status;

  LALMCMCParam *paraHead=NULL;
  REAL4Vector  *step=NULL;
  gsl_matrix *work=NULL; 
  REAL8 aii, aij, ajj;
  INT4 i, j, dim;

  /* set some values */
  dim=parameter->dimension;

  /* draw the mutinormal deviates */
  LALSCreateVector( &status, &step, dim);

  /* copy matrix into workspace and scale it appriopriately */
  work =  gsl_matrix_alloc(dim,dim); 

  gsl_matrix_memcpy( work, covMat );
  gsl_matrix_scale( work, inputMCMC->annealingTemp);
  
  /* check if the matrix if positive definite */
  while ( !XLALCheckPositiveDefinite( work, dim) ) {
    printf("WARNING: Matrix not positive definite!\n");
    /* downweight the off-axis elements */
    for (i=0; i<dim; ++i)
    {
      for (j=0; j<dim; ++j)
      {
        aij=gsl_matrix_get( work, i, j);
        aii=gsl_matrix_get( work, i, i);
        ajj=gsl_matrix_get( work, j, j);  
        
        if ( fabs(aij) > 0.95* sqrt( aii*ajj ) )
        {
          aij=aij/fabs(aij)*0.95*sqrt( aii*ajj );
        }
        gsl_matrix_set( work, i, j, aij);
        gsl_matrix_set( work, j, i, aij);
        printf(" %f", gsl_matrix_get( work, i, j));
      }
      printf("\n");
    }
    exit(0);
  }
    
  /* draw multivariate student distribution with n=2 */
  XLALMultiStudentDeviates( step, work, dim, 2, inputMCMC->randParams); 
  
  /* loop over all parameters */
  for (paraHead=parameter->param,i=0; paraHead; paraHead=paraHead->next,i++)
  { 
  /*  if (inputMCMC->verbose)
      printf("MCMCJUMP: %10s: value: %8.3f  step: %8.3f newVal: %8.3f\n", 
             paraHead->core->name, paraHead->value, step->data[i] , paraHead->value + step->data[i]);*/
    paraHead->value += step->data[i];
	}
  
  XLALMCMCCyclicReflectiveBound(parameter);
  /* destroy the vectors */
  LALSDestroyVector(&status, &step);
  gsl_matrix_free(work);
}

void 
XLALMCMCJumpIntrinsic(
  LALMCMCInput     *inputMCMC,
  LALMCMCParameter *parameter, 
  gsl_matrix       *covMat
  ) 
{ /* </lalVerbatim>  */
  static LALStatus status;

  LALMCMCParam *paraHead=NULL;
  REAL4Vector  *step=NULL;
  gsl_matrix *work=NULL; 
  REAL8 aii, aij, ajj;
  INT4 i, j, dim;

  /* set some values */
  dim=parameter->dimension;

  /* draw the mutinormal deviates */
  LALSCreateVector( &status, &step, dim);
  /* copy matrix into workspace and scale it appriopriately */
  work =  gsl_matrix_alloc(dim,dim); 
  gsl_matrix_memcpy( work, covMat );
  gsl_matrix_scale( work, inputMCMC->annealingTemp);
  
  /* check if the matrix if positive definite */
  while ( !XLALCheckPositiveDefinite( work, dim) ) {
    printf("WARNING: Matrix not positive definite!\n");
    /* downweight the off-axis elements */
    for (i=0; i<dim; ++i)
    {
      for (j=0; j<dim; ++j)
      {
        aij=gsl_matrix_get( work, i, j);
        aii=gsl_matrix_get( work, i, i);
        ajj=gsl_matrix_get( work, j, j);  
        
        if ( fabs(aij) > 0.95* sqrt( aii*ajj ) )
        {
          aij=aij/fabs(aij)*0.95*sqrt( aii*ajj );
        }
        gsl_matrix_set( work, i, j, aij);
        gsl_matrix_set( work, j, i, aij);
        printf(" %f", gsl_matrix_get( work, i, j));
      }
      printf("\n");
    }
    exit(0);
  }

  /* draw multivariate student distribution with n=2 */
  XLALMultiStudentDeviates( step, work, dim, 2, inputMCMC->randParams); 
  
  /* loop over all parameters */
  for (paraHead=parameter->param,i=0; paraHead; paraHead=paraHead->next,i++)
  { 
	if(!strcmp(paraHead->core->name,"long") || !strcmp(paraHead->core->name,"lat")||!strcmp(paraHead->core->name,"time"))
	{;}
  /*  if (inputMCMC->verbose)
      printf("MCMCJUMP: %10s: value: %8.3f  step: %8.3f newVal: %8.3f\n", 
             paraHead->core->name, paraHead->value, step->data[i] , paraHead->value + step->data[i]);*/
    else paraHead->value += step->data[i];
	}
  
  XLALMCMCCyclicReflectiveBound(parameter);
  /* destroy the vectors */
  LALSDestroyVector(&status, &step);
  gsl_matrix_free(work);
}

void XLALMCMCCyclicReflectiveBound(LALMCMCParameter *parameter)
/* Map samples back into parameter space using the simple
cyclic or reflective boundaries - a sampler can use this
function to keep its proposals inside the parameter space */
{
	LALMCMCParam *paraHead=NULL;
	REAL8 delta;
	for (paraHead=parameter->param;paraHead;paraHead=paraHead->next)
	{
		if(paraHead->core->wrapping) /* For cyclic boundaries */
		{
			delta = paraHead->core->maxVal - paraHead->core->minVal;
			while ( paraHead->value > paraHead->core->maxVal) 
				paraHead->value -= delta;
			while ( paraHead->value < paraHead->core->minVal) 
			paraHead->value += delta;
		}
		else /* Use reflective boundaries */
		{
			if(paraHead->core->maxVal < paraHead->value) paraHead->value-=2.0*(paraHead->value - paraHead->core->maxVal);
			if(paraHead->core->minVal > paraHead->value) paraHead->value+=2.0*(paraHead->core->minVal - paraHead->value);
		}
	}
}

INT4 XLALMCMCCheckParameter(
			   LALMCMCParameter *parameter,
			   const char *name)
{
  /* Check for existance of name in parameter */
  LALMCMCParam *param=NULL;
  param=parameter->param;
  while(param) {if(!strcmp(param->core->name, name)) return 1; else param=param->next;}
  return 0;

}


/* *****************************
XLALMCMCAddParam
  ***************************** */
/*  <lalVerbatim file="XLALMCMCAddParamCP"> */
void
XLALMCMCAddParam(
  LALMCMCParameter  *parameter,
  const char        *name, 
  REAL8              value,
  REAL8              minValue,
  REAL8              maxValue,
  INT4               wrapping
  )
{ /* </lalVerbatim>  */

  LALMCMCParam* paraPointer;

  if ( !parameter )
  {
    fprintf( stderr, "ERROR in XLALMCMCAddParam: 'parameter' is a NULL pointer\n");
    exit(0);
  }

  if ( !parameter->param )
  {
    paraPointer = parameter->param = (LALMCMCParam*) LALMalloc( sizeof(LALMCMCParam) );
  }
  else
  {

    /* first search the end of the line */  
    paraPointer=parameter->param;
    while ( paraPointer->next )
    {
      paraPointer=paraPointer->next;
    }
    
    paraPointer = paraPointer->next = (LALMCMCParam*) LALMalloc( sizeof(LALMCMCParam) );
  }

  /* set the next pointer to NULL always */
  paraPointer->next = NULL;
  
  /* allocate the sub structure (once!!) */
  paraPointer->core = (LALMCMCSubParam*) LALMalloc( sizeof(LALMCMCSubParam) );

  /* fill the sub-structure */
  memcpy( paraPointer->core->name, name, 30 );
  paraPointer->core->minVal = minValue;
  paraPointer->core->maxVal = maxValue;
  paraPointer->core->wrapping = wrapping;
  paraPointer->core->chain = NULL;
  paraPointer->value = value;

  /*printf("MCMCInit parameter %s with range %f - %f - %f\n", name, minValue, value, maxValue ); */

  /* increase dimension */
  parameter->dimension++;  
}



/* *****************************
XLALMCMCGetParam
  ***************************** */
/*  <lalVerbatim file="XLALMCMCGetParamCP"> */
LALMCMCParam*
XLALMCMCGetParam(
  LALMCMCParameter* parameter,
  const char* name
  )
{ /* </lalVerbatim>  */

  LALMCMCParam* param=NULL;

  param=parameter->param;
  while (param)
  {
    /* check if the name is correct */
    if (!strcmp(param->core->name, name))
    {
      return param;
    }

    /* go to the next instance */
    param=param->next;
  }

  /* parameter not found .... */
  fprintf( stderr, 
           "WARNING: parameter '%s' unknown!\n",name );

  return param;
}


/* *****************************
XLALMCMCGetParameter
  ***************************** */
/*  <lalVerbatim file="XLALMCMCGetParameterCP"> */
REAL8
XLALMCMCGetParameter(
  LALMCMCParameter* parameter,
  const char* name
  )
{ /* </lalVerbatim>  */

  LALMCMCParam* param=NULL;

  param=XLALMCMCGetParam( parameter, name);
  
  if (param)
  {
    return param->value;
  }

  return 0.0;
}


/* *****************************
XLALMCMCSetParameter
  ***************************** */
/*  <lalVerbatim file="XLALMCMCSetParameterCP"> */
void
XLALMCMCSetParameter(
  LALMCMCParameter* parameter,
  const char* name,
  REAL8 value
  )
{ /* </lalVerbatim>  */

  LALMCMCParam* param;
  param=XLALMCMCGetParam( parameter, name);

  if (param)
  {
    param->value=value;
  }
  
}


/* *****************************
XLALMCMCCopyPara
  ***************************** */
/*  <lalVerbatim file="XLALMCMCCopyParaCP"> */
void
XLALMCMCCopyPara(
  LALMCMCParameter **parameterOutPtr,
  LALMCMCParameter *parameterIn
  )
{ /* </lalVerbatim>  */

  /* deep copy of the param structure, parameterOut must be a already allocated pointer */
  LALMCMCParam *outPointer= NULL;
  LALMCMCParam *inPointer = NULL;
  LALMCMCParam *before    = NULL;
  
  LALMCMCParameter *parameterOut = NULL;

  parameterOut = *parameterOutPtr;

  /* check if pointer is allocated */
  if ( !parameterOut )
  {
    fprintf(stderr," ERROR in XLALMCMCCopyPara: 'parameter' not allocated");
    exit(0);
  }
  
  /* deep copy of the upper structure */
  /*memcpy( parameterOut, parameterIn, sizeof(LALMCMCParameter) );*/
  parameterOut->dimension     = parameterIn->dimension;
  parameterOut->logLikelihood = parameterIn->logLikelihood;
  parameterOut->logPrior      = parameterIn->logPrior;

   
  /* allocate the first of the param-pointers if needed */
  if(parameterOut->param==NULL) parameterOut->param = (LALMCMCParam*) LALMalloc( sizeof(LALMCMCParam) );
 
  outPointer = parameterOut->param;
  /*outPointer->next = NULL; */

  /* loop over all the param pointers in the in-(mother) structure */
  for ( inPointer=parameterIn->param; inPointer; inPointer = inPointer->next )
  {
    /* check if the next-pointer is already allocated (happens only the first time) */
    if ( !outPointer ) 
    {
      outPointer = (LALMCMCParam*) LALMalloc( sizeof(LALMCMCParam) );
      before->next = outPointer;           
      outPointer->next = NULL; 
    }

    /* copy the sub-param structure: just copy the pointer ... */
    outPointer->core = inPointer->core;

    /* ... and just the actual value corresponding to this parameter */
    outPointer->value = inPointer->value;

    before = outPointer;      /* remember the previous cell */
    outPointer = outPointer->next;
  }


}


/* *****************************
XLALMCMCFreePara
  ***************************** */
/*  <lalVerbatim file="XLALMCMCFreeParaCP"> */
void
XLALMCMCFreePara(
  LALMCMCParameter *parameter
  )
{ /* </lalVerbatim>  */

  LALMCMCParam* param=NULL; 
  LALMCMCParam* thisParam=NULL;

  
  param=parameter->param;
  while (param)
  {
    thisParam = param;
    param = param->next;
	LALFree( thisParam->core);
    LALFree( thisParam );    
  }
  parameter->param = NULL;
  parameter->dimension=0;
}


/* *****************************
XLALMCMCDestroyPara
  ***************************** */
/*  <lalVerbatim file="XLALMCMCDestroyParaCP"> */
void
XLALMCMCDestroyPara(
  LALMCMCParameter **parameter
  )
{ /* </lalVerbatim>  */

  LALMCMCParameter* para;
  LALMCMCParam* param=NULL; 
  LALMCMCParam* paramNext=NULL;
  
  para=*parameter;

  param=para->param;
  while (param)
  {
    paramNext=param->next;
    LALFree( param->core );
    LALFree( param );
    param=paramNext;
    
  }

  LALFree( para );
}


/* *****************************
XLALMultiStudentDeviates
  ***************************** */
/*  <lalVerbatim file="XLALMultiStudentDeviatesCP"> */
void
XLALMultiStudentDeviates( 
  REAL4Vector  *vector,
  gsl_matrix   *matrix,
  UINT4         dim,
  UINT4         n,
  RandomParams *randParam
  )
{ /* </lalVerbatim> */
  static const char *func = "LALMultiStudentDeviates";

  static LALStatus status;

  REAL4Vector *dummy=NULL;
  REAL4 chi=0.0, factor;
  UINT4 i;

  /* check input arguments */
  if (!vector || !matrix || !randParam)
    XLAL_ERROR_VOID( func, XLAL_EFAULT );
  
  if (dim<1)
    XLAL_ERROR_VOID( func, XLAL_EINVAL );

  if (n<1)
    XLAL_ERROR_VOID( func, XLAL_EINVAL );


   /* first draw from MVN */
  XLALMultiNormalDeviates( vector, matrix, dim, randParam);


  /* then draw from chi-square with n degrees of freedom;
     this is the sum d_i*d_i with d_i drawn from a normal 
     distribution. */
  LALSCreateVector( &status, &dummy, n);
  LALNormalDeviates( &status, dummy, randParam);

  /* calculate the chisquare distributed value */
  for (i=0; i<n; i++) 
  {
    chi+=dummy->data[i]*dummy->data[i];
  }

  /* destroy the helping vector */
  LALSDestroyVector( &status, &dummy );

  /* now, finally, calculate the distribution value */
  factor=sqrt(n/chi);
  for (i=0; i<dim; i++) 
  {
    vector->data[i]*=factor;
  }

}


/* Reference: http://www.mail-archive.com/help-gsl@gnu.org/msg00631.html*/
/*  <lalVerbatim file="XLALMultiNormalDeviatesCP"> */
void
XLALMultiNormalDeviates( 
  REAL4Vector *vector, 
	gsl_matrix *matrix, 
  UINT4 dim, 
  RandomParams *randParam
  )
{/* </lalVerbatim> */
  static LALStatus status;

  UINT4 i=0;
  gsl_matrix *work=NULL;
  gsl_vector *result = NULL;
  
  static const char *func = "LALMultiNormalDeviates";
  
  /* check input arguments */
  if (!vector || !matrix || !randParam)
    XLAL_ERROR_VOID( func, XLAL_EFAULT );
  
  if (dim<1)
    XLAL_ERROR_VOID( func, XLAL_EINVAL );

  /* copy matrix into workspace */
  work =  gsl_matrix_alloc(dim,dim); 
  gsl_matrix_memcpy( work, matrix );

  /* compute the cholesky decomposition */
  gsl_linalg_cholesky_decomp(work);
  
  /* retrieve the normal distributed random numbers (LAL procedure) */
  LALNormalDeviates( &status, vector, randParam);

  /* store this into a gsl vector */
  result = gsl_vector_alloc ( (int)dim );
  for (i = 0; i < dim; i++)
  {
    gsl_vector_set (result, i, vector->data[i]);
  }

  /* compute the matrix-vector multiplication */
  gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);

  /* recopy the results */
  for (i = 0; i < dim; i++)
  {
    vector->data[i]=gsl_vector_get (result, i);
  }

  /* free unused stuff */
  gsl_matrix_free(work);
  gsl_vector_free(result);
  
}
 
 
/*  <lalVerbatim file="XLALCheckPositiveDefiniteCP"> */
UINT4
XLALCheckPositiveDefinite( 
  gsl_matrix       *matrix,
  UINT4            dim
  )
{/* </lalVerbatim> */
  gsl_matrix  *m     = NULL;
  gsl_vector  *eigen = NULL;
  gsl_eigen_symm_workspace *workspace = NULL;
  UINT4 i;
  
  /* copy input matrix */
  m =  gsl_matrix_alloc( dim,dim ); 
  gsl_matrix_memcpy( m, matrix);  
  
  /* prepare variables */
  eigen = gsl_vector_alloc ( dim );
  workspace = gsl_eigen_symm_alloc ( dim );
  
  /* compute the eigen values */
  gsl_eigen_symm ( m,  eigen, workspace );
  
  /* test the result */
  for (i = 0; i < dim; i++)
    {
      /* printf("diag: %f | eigen[%d]= %f\n", gsl_matrix_get( matrix,i,i), i, eigen->data[i]);*/
    if (eigen->data[i]<0) 
    {
      printf("NEGATIVE EIGEN VALUE!!! PANIC\n");
      return 0;
    }
  }

  /* freeing unused stuff */
  gsl_eigen_symm_free( workspace);
  gsl_matrix_free(m);
  gsl_vector_free(eigen);
  
  return 1;
}
 
