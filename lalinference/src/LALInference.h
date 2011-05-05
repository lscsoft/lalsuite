/*
 *
 *  LALInference:             Bayesian Followup        
 *  include/LALInference.h:   main header file
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch
 *
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
 * \file LALInference.h
 * \brief Main header file for LALInference common routines
 */

#ifndef LALInference_h
#define LALInference_h

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

#define VARNAME_MAX 128

# include <lal/LALStdlib.h>
# include <lal/LALConstants.h>
# include <lal/SimulateCoherentGW.h>
# include <lal/GeneratePPNInspiral.h>
# include <lal/LIGOMetadataTables.h>
# include <lal/LALDatatypes.h>
# include <lal/FindChirp.h>
# include <lal/Window.h>
#include <lal/LALString.h>

#include <lal/SFTutils.h>
#include <lal/SFTfileIO.h>
#include <lal/LALDetectors.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/BinaryPulsarTiming.h>


#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>

//...other includes

struct tagLALInferenceRunState;
struct tagLALInferenceIFOData;

/*Data storage type definitions*/

/** An enumerated type for denoting the type of a variable
*/
typedef enum {
  INT4_t, 
  INT8_t,
  UINT4_t,
  REAL4_t, 
  REAL8_t, 
  COMPLEX8_t, 
  COMPLEX16_t, 
  gslMatrix_t,
  REAL8Vector_t,
  UINT4Vector_t,
  string_t
} LALInferenceVariableType;

/** An enumerated type for denoting the domain
*/
typedef enum {
  timeDomain, 
  frequencyDomain
} LALInferenceDomain;

/** An enumerated type for denoting the topolology of a parameter.
*/
typedef enum {
	PARAM_LINEAR, /** A parameter that simply has a maximum and a minimum */
	PARAM_CIRCULAR, /** A parameter that is cyclic, such as an angle between 0 and 2\pi */
	PARAM_FIXED,    /** A parameter that never changes, functions should respect this */
	PARAM_OUTPUT    /** A parameter changed by an inner code and passed out */
} LALInferenceParamVaryType;

/** An enumerated type for denoting a type of taper
*/
typedef enum
{
	INFERENCE_TAPER_NONE,
	INFERENCE_TAPER_START,
	INFERENCE_TAPER_END,
	INFERENCE_TAPER_STARTEND,
	INFERENCE_TAPER_NUM_OPTS,
	INFERENCE_RING,
	INFERENCE_SMOOTH
}  LALInferenceApplyTaper;


extern size_t typeSize[];

/** The LALInferenceVariableItem list node structure
 * This should only be accessed using the accessor functions below
 * Implementation may change to hash table
*/
typedef struct
tagVariableItem
{
  char                    name[VARNAME_MAX];
  void                    *value;
  LALInferenceVariableType            type;
  LALInferenceParamVaryType			  vary;
  struct tagVariableItem  *next;
} LALInferenceVariableItem;

/** The LALInferenceVariables list structure
 * Should only be accessed using the accessor functions below
 */
typedef struct
tagLALInferenceVariables
{
  LALInferenceVariableItem * head;
  INT4 dimension;
} LALInferenceVariables;


const char *translateInternalToExternalParamName(const char *inName);
int fprintParameterNonFixedHeaders(FILE *out, LALInferenceVariables *params);

/** Return a pointer to the variable asked for */
void *LALInferenceGetVariable(LALInferenceVariables * vars, const char * name);
/** Get number of dimensions in this variable */
INT4 LALInferenceGetVariableDimension(LALInferenceVariables *vars);
/** Get number of dimensions which are not fixed to a certain value */
INT4 LALInferenceGetVariableDimensionNonFixed(LALInferenceVariables *vars);

 
LALInferenceVariableType LALInferenceGetVariableTypeByIndex(LALInferenceVariables *vars, int index);


/** Get the type of the variable, returns LALInferenceVariableType (see above) */
LALInferenceVariableType LALInferenceGetVariableType(LALInferenceVariables *vars, const char *name);

/** Get the LALInferenceParamVaryType of the variable, see types above */
LALInferenceParamVaryType LALInferenceGetVariableVaryType(LALInferenceVariables *vars, const char *name);

/** Get the name of the index-th variable */
char *LALInferenceGetVariableName(LALInferenceVariables *vars, int index);

/** Set a variable with a value, pass a void * to the value as argument */
void LALInferenceSetVariable(LALInferenceVariables * vars, const char * name, void * value);

/** Add a variable to the list */
void LALInferenceAddVariable(LALInferenceVariables * vars, const char * name, void * value, 
	LALInferenceVariableType type, LALInferenceParamVaryType vary);

/** Remove a variable fm the list */
void LALInferenceRemoveVariable(LALInferenceVariables *vars,const char *name);

/** Checks for a variable being present in the list, returns 1(==true) or 0 */
int  LALInferenceCheckVariable(LALInferenceVariables *vars,const char *name);

/** Delete the variables in this structure (does not free the LALInferenceVariables itself) */
void LALInferenceDestroyVariables(LALInferenceVariables *vars);

/** Deep copy the variables from one to another LALInferenceVariables structure */
void LALInferenceCopyVariables(LALInferenceVariables *origin, LALInferenceVariables *target);

/** Print variables to stdout */
void LALInferencePrintVariables(LALInferenceVariables *var);

/** Check for equality in two variables */
int LALInferenceCompareVariables(LALInferenceVariables *var1, LALInferenceVariables *var2);


//Wrapper for template computation 
//(relies on LAL libraries for implementation) <- could be a #DEFINE ?
//typedef void (LALTemplateFunction) (LALInferenceVariables *currentParams, struct tagLALInferenceIFOData *data); //Parameter Set is modelParams of LALInferenceIFOData
/** Type declaration for template function */
typedef void (LALInferenceTemplateFunction) (struct tagLALInferenceIFOData *data);


//Jump proposal distribution
//Computes proposedParams based on currentParams and additional variables
//stored as proposalArgs, which could include correlation matrix, etc.,
//as well as forward and reverse proposal probability.
//A jump proposal distribution function could call other jump proposal
//distribution functions with various probabilities to allow for multiple
//jump proposal distributions
/** Type declaration for Proposal function */
typedef void (LALInferenceProposalFunction) (struct tagLALInferenceRunState *runState,
	LALInferenceVariables *proposedParams);

/** Type declaration for prior function */
typedef REAL8 (LALInferencePriorFunction) (struct tagLALInferenceRunState *runState,
	LALInferenceVariables *params);

//Likelihood calculator 
//Should take care to perform expensive evaluation of h+ and hx 
//only once if possible, unless necessary because different IFOs 
//have different data lengths or sampling rates 
/** Type declaration for likelihood function */
typedef REAL8 (LALInferenceLikelihoodFunction) (LALInferenceVariables *currentParams,
        struct tagLALInferenceIFOData * data, LALInferenceTemplateFunction *template);

/** Perform one step of an algorithm, replaces runState->currentParams */
typedef void (LALInferenceEvolveOneStepFunction) (struct tagLALInferenceRunState *runState);

/** Type declaration for an algorithm function; will distinguish MCMC from NestedSampling, etc */
typedef void (LALInferenceAlgorithm) (struct tagLALInferenceRunState *runState);

/** Structure containing inference run state
 * This includes pointers to the function types required to run
 * the algorithm, and data structures as required */
typedef struct 
tagLALInferenceRunState
{
  ProcessParamsTable        *commandLine; /** A ProcessParamsTable with command line arguments */
  LALInferenceAlgorithm              *algorithm; /** The algorithm function */
  LALInferenceEvolveOneStepFunction  *evolve; /** The algorithm's single iteration function */
  LALInferencePriorFunction          *prior; /** The prior for the parameters */
  LALInferenceLikelihoodFunction     *likelihood; /** The likelihood function */
  LALInferenceProposalFunction       *proposal; /** The proposal function */
  LALInferenceTemplateFunction       *template; /** The template generation function */
  struct tagLALInferenceIFOData      *data; /** The data from the interferometers */
  LALInferenceVariables              *currentParams, /** The current parameters */
    *priorArgs,                                      /** Any special arguments for the prior function */
    *proposalArgs,                                   /** Any special arguments for the proposal function */
    *algorithmParams;                                /** Parameters which control the running of the algorithm*/
  LALInferenceVariables				**livePoints; /** Array of live points for Nested Sampling */
  LALInferenceVariables **differentialPoints;        /** Array of points for differential evolution */
  size_t differentialPointsLength;                   /** This should be removed can be given as an algorithmParams entry */
  REAL8			currentLikelihood;  /** This should be removed, can be given as an algorithmParams or proposalParams entry */
  REAL8                 currentPrior;       /** This should be removed, can be given as an algorithmParams entry */
  gsl_rng               *GSLrandom;         /** A pointer to a GSL random number generator */
} LALInferenceRunState;


/* LALInferenceRunState *initialize(ProcessParamsTable *commandLine); */

#define DETNAMELEN 256

/** Structure to contain IFO data.
 *  Some fields may be left empty if not needed
*/
typedef struct
tagLALInferenceIFOData
{
  char                       name[DETNAMELEN]; /** Detector name */
  REAL8TimeSeries           *timeData,         /** A time series from the detector */
                            *timeModelhPlus, *timeModelhCross, /** Time series model buffers */
                            *whiteTimeData, *windowedTimeData; /** white is not really white, but over-white. */
  /* Stores the log(L) for the model in presence of data.  These were
     added to allow for individual-detector log(L) output.  The
     convention is that loglikelihood always stores the log(L) for the
     model in freqModel... or timeModel....  When a jump is accepted,
     that value is copied into acceptedloglikelihood, which is the
     quantity that is actually output in the output files. */
  REAL8                      nullloglikelihood, loglikelihood, acceptedloglikelihood; 
  REAL8                      fPlus, fCross; /** Detector responses */
  REAL8                      timeshift;     /** What is this? */
  COMPLEX16FrequencySeries  *freqData,      /** Buffer for frequency domain data */
                            *freqModelhPlus, *freqModelhCross, /** Buffers for frequency domain models */
                            *whiteFreqData; /* Over-white. */
  COMPLEX16TimeSeries       *compTimeData, *compModelData; /** Complex time series data buffers */
  LIGOTimeGPSVector         *dataTimes;                    /** Vector of time stamps for time domain data */
  LALInferenceVariables              *modelParams;         /** Parameters used when filling the buffers - template functions should copy to here */
  LALInferenceVariables		    *dataParams; /* Optional data parameters */
  LALInferenceDomain                 modelDomain;         /** Domain of model */
  REAL8FrequencySeries      *oneSidedNoisePowerSpectrum;  /** one-sided Noise Power Spectrum */
  REAL8TimeSeries           *timeDomainNoiseWeights; /** Roughly, InvFFT(1/Noise PSD). */
  REAL8Window               *window;                 /** A window */
  REAL8FFTPlan              *timeToFreqFFTPlan, *freqToTimeFFTPlan; /** Pre-calculated FFT plans for forward and reverse FFTs */
  REAL8                     fLow, fHigh;	/** integration limits for overlap integral in F-domain */
  LALDetector               *detector;          /** LALDetector structure for where this data came from */
  BarycenterInput           *bary;              /** Barycenter information */
  EphemerisData             *ephem;             /** Ephemeris data */
  LIGOTimeGPS		    epoch;              /** The epoch of this observation (the time of the first sample) */

  struct tagLALInferenceIFOData      *next;     /** A pointer to the next set of data for linked list */
} LALInferenceIFOData;

/** Returns the element of the process params table with "name" */
ProcessParamsTable *LALInferenceGetProcParamVal(ProcessParamsTable *procparams,const char *name);

/** What is this? */
void LALInferenceParseCharacterOptionString(char *input, char **strings[], UINT4 *n);

/** Return a ProcessParamsTable from the command line arguments */
ProcessParamsTable *LALInferenceParseCommandLine(int argc, char *argv[]);

/** Output the command line based on the ProcessParamsTable */
void LALInferencePrintCommandLine(ProcessParamsTable *procparams, char *str);

/** This should be moved elsewhere */
void COMPLEX16VectorSubtract(COMPLEX16Vector * out, const COMPLEX16Vector * in1, const COMPLEX16Vector * in2);

/** Execute FFT and inverse FFT for the data */
void LALInferenceExecuteFT(LALInferenceIFOData *IFOdata);
void LALInferenceExecuteInvFT(LALInferenceIFOData *IFOdata);

/** This should be removed */
void die(const char *message);

/** Return the list node for "name" - do not rely on this */
LALInferenceVariableItem *LALInferenceGetItem(LALInferenceVariables *vars,const char *name);
/** Return the list node for the index-th item - do not rely on this */
LALInferenceVariableItem *LALInferenceGetItemNr(LALInferenceVariables *vars, int index);

/** Output the sample to file *fp, in ASCII format */
void LALInferencePrintSample(FILE *fp,LALInferenceVariables *sample);
/** Output only non-fixed parameters */
void LALInferencePrintSampleNonFixed(FILE *fp,LALInferenceVariables *sample);

void mc2masses(double mc, double eta, double *m1, double *m2);

/* Differential Evolution and Common Format Posterior File Parsing
   Utilities. */

/** Returns an array of header strings (terminated by NULL) from a common-format output file */
char **getHeaderLine(FILE *inp);

/** Turns common-format column names into our internal parameter
   names. */
char *colNameToParamName(const char *colName);

/** Reads one line from the given file and stores the values there into
   the variable structure, using the given header array to name the
   columns.  Returns 0 on success. */
int processParamLine(FILE *inp, char **headers, LALInferenceVariables *vars);

#endif

