/*******************************************************************************
  Matt Pitkin - 08/08/07

  pulsar_parameter_estimation.h

  Header file for pulsar_parameter_estimation.c

*******************************************************************************/

/*
  Author:
*/

#ifndef _PULSAR_PARAMETER_ESTIMATION_H
#define _PULSAR_PARAMETER_ESTIMATION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/time.h>

/* LAL headers */
#include <lal/LALStdlib.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALDatatypes.h>
#include <lal/PulsarDataTypes.h>
#include <lal/AVFactories.h>
#include <lal/SkyCoordinates.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/BinaryPulsarTiming.h>
#include <lal/Random.h>
#include <lal/LALString.h>
#include <lal/SFTutils.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/MatrixUtils.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/CWFastMath.h>
#include <lal/TimeSeries.h>

#include <lalapps.h>

#ifdef __cplusplus
extern "C" {
#endif

/** define macros */
#define ROUND(x) (floor(x+0.5))

/* macro to perform logarithmic addition log(exp(x) + exp(y)) */
#define PLUS(x,y) ( x>y ? x+log(1.+exp(y-x)) : y+log(1.+exp(x-y)) )

/* macro to calculate the area of a trapezium for on logs */
/* logarea = log(0.5) + log(width) + log(exp(logHeight1) + exp(logHeight2)) */
#define LOG_TRAPEZIUM(h1, h2, w) ( -0.693147180559945 + log(w) + PLUS(h1, h2) )

#define MAXLENGTH 1000000
#define MAXPARAMS 35

#define SIXTH 0.166666666666666666666666666666666666666667L
#define TWENTYFOURTH 0.04166666666666666666666666666666666666666667L

/** define structures */

/* a structure to hold the four instrinsic pulsar parameters */
typedef struct tagIntrinsicPulsarVariables{
  REAL8 h0;         /* GW amplitude */
  REAL8 phi0;       /* initial GW phase */
  REAL8 ci;         /* cos(inclination angle) */
  REAL8 psi;        /* polarisation angle */

  /* derived variables to allow quicker computation on the grid */
  /* MUST BE DEFINED PRIOR TO THE LIKELIHOOD CALCULATION */
  REAL8 Xplus;      /* 0.5*(1+cos^2(iota)) */
  REAL8 Xcross;     /* cos(iota) */
  REAL8 Xpsinphi_2; /* 0.5*Xplus*sin(phi) */
  REAL8 Xcsinphi_2; /* 0.5*Xcross*sin(phi) */
  REAL8 Xpcosphi_2; /* 0.5*Xplus*cos(phi) */
  REAL8 Xccosphi_2; /* 0.5*Xcross*cos(phi) */

  REAL8 Aplus;      /* 0.5*h0*(1+cos^2(iota)) */
  REAL8 Across;     /* h0*cos(iota) */
  REAL8 Apsinphi_2; /* 0.5*Aplus*sin(phi) */
  REAL8 Acsinphi_2; /* 0.5*Across*sin(phi) */
  REAL8 Apcosphi_2; /* 0.5*Aplus*cos(phi) */
  REAL8 Accosphi_2; /* 0.5*Across*cos(phi) */
}IntrinsicPulsarVariables;

typedef struct tagMeshGrid{
  IntrinsicPulsarVariables minVals;
  IntrinsicPulsarVariables maxVals;
  IntrinsicPulsarVariables delta; /* parameter step sizes */

  INT4 h0Steps, phiSteps, psiSteps, ciotaSteps;

  /* the number of steps in the range for the lookup table */
  INT4 psiRangeSteps;
  INT4 timeRangeSteps;
}MeshGrid;

/* a structure to hold the MCMC initial parameters */
typedef struct tagMCMCParams{
  /* MCMC parameters */
  INT4 doMCMC;                     /* flag can be true 1 or false 0 */
  INT4 iterations;                 /* length of MCMC chain */
  INT4 burnIn;                     /* length of burn in period */
  INT4 outputBI;                   /* output the burn in chain */
  REAL8 temperature;               /* temperature of simulated annealing during
                                      burn in */
  IntrinsicPulsarVariables sigmas; /* standard deviations of the proposal
                                      distributions */
  REAL8 h0scale;                   /* scale factor to multiply h0 sigma by */

  INT4 outputRate;                 /* rate at which to output chain e.g. 10
                                      means output every tenth sample */

  CHAR glitchTimes[500];           /* string containing the times of glitches
                                      in MJD */
  INT4 nGlitches;                  /* the number of glitches */
  REAL8 glitchCut;                 /* the number of seconds of data to ignore
                                      around the time of the glitch */
}MCMCParams;

typedef struct tagPriorVals{
  CHAR *h0Prior, *phiPrior, *psiPrior, *iotaPrior; /* prior type e.g.
"uniform", "jeffreys" or "gaussian" */
  CHAR *priorFile; /* a file containing the h0xci prior distribution */
  REAL8Vector *h0vals, *civals;
  REAL8 **h0cipdf;
  REAL8 minh0ci, maxh0ci; /* minimum and maximum prior values */
  
  IntrinsicPulsarVariables vars;

  /* for gaussian prior pass in the mean and standard deviation */
  REAL8 meanh0, stdh0;
  REAL8 meanphi, stdphi;
  REAL8 meanpsi, stdpsi;
  REAL8 meaniota, stdiota;
}PriorVals;

/* a structure to hold the command line input values */
typedef struct tagInputParams{
  CHAR detectors[40];
  CHAR pulsar[20];
  CHAR parFile[256];   /* pulsar parameter file */
  CHAR *matrixFile;    /* pulsar parameter covariance matrix file */
  UINT4 usecov;        /* set whether or not to use a covariance matrix prior */
  
  CHAR inputDir[256];
  CHAR outputDir[256];

  LALSource psr; /* pulsar position values */

  /* prior values */
  INT4 usepriors;
  PriorVals priors;

  /* grid values */
  MeshGrid mesh;

  /* maximum and minumum length of time over which we assume stationarity of
     the data e.g. min 5 mins and max 30 mins */
  INT4 chunkMin;
  INT4 chunkMax;

  /* parameters for an MCMC */
  MCMCParams mcmc;

  REAL8 dob;                      /* degree of belief for upper limit */

  INT4 outputPost; /* flag for whether of not to output the full posterior */

  CHAR earthfile[256];
  CHAR sunfile[256];
  CHAR timefile[256];

  INT4 onlyjoint;

  INT4 verbose;
}InputParams;

/* a detector response function lookup table structure */
typedef struct tagDetRespLookupTable{
  LALDetAMResponse **lookupTable; /* array containing the lookup table */
  INT4 psiSteps;                   /* number of steps across psi range */
  INT4 timeSteps;                  /* number of steps across a day */
}DetRespLookupTable;

typedef struct tagDataStructure{
  COMPLEX16Vector *data;     /* real and imaginary B_ks */
  REAL8Vector *times;        /* time stamp for each B_k */
  INT4Vector *chunkLengths;  /* length of each data chunk e.g. 30 mins */
  INT4 chunkMax;              /* maximum chunk length to be used e.g. 30 mins */
  INT4 chunkMin;              /* minumum chunk length to be used e.g. 5 mins */

  REAL8Vector *sumData;      /* sum over the data for each chunk (Re(B)^2 +
                               Im(B)^2)*/
  REAL8Vector *sumModel;     /* sum over the model for each chunk (Re(y)^2 +
                               Im(y)^2) not including the h0 scaling */
  REAL8Vector *sumDataModel; /* sum over the model and data for each chunk
                               (Re(B)*Re(y) + Im(B)*Im(y)) */
  DetRespLookupTable *lookupTable;
}DataStructure;

typedef struct tagOutputParams{
  CHAR *outputDir;  /* directory to output files to */
  CHAR *psr;        /* pulsar name */
  CHAR *det;        /* detector e.g. H1, H2, L1, G1, V1, Joint */

  CHAR *margParam;  /* parameter over which to not marginalise */

  INT4 outPost;      /* flag for whether to output the full posterior 1=yes,
                       0=no */
  REAL8 dob;       /* the degree-of-belief for an upper limits calculation -
                       don't perform calc if set to zero */ 
}OutputParams;


typedef struct tagResults{
  REAL8 evidence;        /* the evidence */
  REAL8 h0UpperLimit;    /* the h0 upper limit */
  REAL8 ellipticity;     /* the ellipticity upper limit */
  REAL8 spinDown;        /* the ratio to the spin-down UL */
}Results;

typedef struct tagParamData{
  const CHAR *name;  /* parameter name as given by the conventions in the .mat file
(see param variable in TEMPOs mxprt.f file */
  REAL8 val;   /* parameter value */
  REAL8 sigma; /* standard deviation on the parameter as read from the .par
file */
  INT4 matPos;  /* row position in the covariance matrix */
}ParamData;

/** define functions */

/* function to get the input arguments from the command line */
void get_input_args(InputParams *inputParams, INT4 argc, CHAR *argv[]);

/* function to allocate memory for a likelihood array - it will be index as
   logLike[phi][cosiota][psi][h0] */
REAL8 ****allocate_likelihood_memory(MeshGrid mesh);

/* function to create a log likelihood array over the parameter grid */
REAL8 create_likelihood_grid(DataStructure data, REAL8 ****logLike, 
  MeshGrid mesh);

/* a function to compute the log likelihood of the data given some parameters 
- this function loops over the h0 part of the likelihood array for speed */
REAL8 log_likelihood(REAL8 *likeArray, DataStructure data,
  IntrinsicPulsarVariables vars, MeshGrid mesh, REAL8Vector *dphi);

/* a function to sum over the data */
void sum_data(DataStructure data);

/* a function to give the log prior given a set of prior ranges and parameter
values */
REAL8 log_prior(PriorVals priors, MeshGrid mesh);

/* function to compute the log posterior */
REAL8 log_posterior(REAL8 ****logLike, PriorVals prior, MeshGrid mesh,
  OutputParams output);

/* marginalise posterior over requested parameter and output the Results if
requested */
Results marginalise_posterior(REAL8 ****logPost, MeshGrid mesh, 
  OutputParams output);

/* detector response lookup table function  - this function will output a lookup
table of points in time and psi, covering a day from the start time (t0) and
from -pi/4 to pi/4 in psi */
void response_lookup_table(REAL8 t0, LALDetAndSource detAndSource,
  DetRespLookupTable *lookupTable);

/* factorial function */
REAL8 log_factorial(INT4 num);

/* function to combine log likelihoods to give a joint likelihood */
void combine_likelihoods(REAL8 ****logLike1, REAL8 ****logLike2, 
  MeshGrid mesh);

/* function to calculate the upper limit - use quadratic spline interpolation 
  between points around the upper limit */
REAL8 get_upper_limit(REAL8 *cumsum, REAL8 limit, MeshGrid mesh);

/* function to perform the MCMC parameter estimation */
void perform_mcmc(DataStructure *data, InputParams input, INT4 numDets, 
  CHAR *det, LALDetector *detpos, EphemerisData *edat,
  TimeCorrectionData *tdat, TimeCorrectionType ttype);

/* function to get a vector of phases for all the points in data */
REAL8Vector *get_phi(DataStructure data, BinaryPulsarParams params,
  BarycenterInput bary, EphemerisData *edat, TimeCorrectionData *tdat,
  TimeCorrectionType ttype);

/* function to get the lengths of consecutive chunks of data */
void get_chunk_lengths(DataStructure data);

REAL8Array *cholesky_decomp( REAL8Array *M, const CHAR* uOrl );

REAL8Array *read_correlation_matrix( CHAR *matrixFile, 
  BinaryPulsarParams params, ParamData *data );

REAL8Array *create_covariance_matrix( ParamData *data, REAL8Array *corMat, 
  INT4 isinv );

REAL8Array *check_positive_definite( REAL8Array *matrix );

REAL8Array *convert_to_positive_definite( REAL8Array *nonposdef );

ParamData *multivariate_normal_deviates( REAL8Array *covmat, ParamData *means,
  RandomParams *randomParams );

void set_mcmc_pulsar_params( BinaryPulsarParams *pulsarParams, ParamData *data);

#ifdef __cplusplus
}
#endif

#endif /* _PULSAR_PARAMETER_ESTIMATION_H */
