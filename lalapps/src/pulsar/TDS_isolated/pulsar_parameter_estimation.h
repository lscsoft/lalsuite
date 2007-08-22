/*******************************************************************************
  Matt Pitkin - 08/08/07
  
  pulsar_parameter_estimation.h
  
  Header file for pulsar_parameter_estimation.c
  
*******************************************************************************/

/* 
  Author:
  $Id$
*/

#ifndef _PULSAR_PARAMETER_ESTIMATION_H
#define _PULSAR_PARAMETER_ESTIMATION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>

/* LAL headers */
#include <lal/LALStdlib.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/SkyCoordinates.h> 
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/BinaryPulsarTiming.h>

#ifdef  __cplusplus
extern "C" {
#endif


/** define macros */
#define ROUND(x) (floor(x+0.5))

/* macro to perform logarithmic addition log(exp(x) + exp(y)) */
#define PLUS(x,y) ( x>y ? x+log(1.+exp(y-x)) : y+log(1.+exp(x-y)) )

/* Usage format string */
#define USAGE \
"Usage: %s [options]\n\n"\
" --help              display this message\n"\
" --verbose           display all error messages\n"\
" --detectors         all IFOs with data to be analysed e.g. \"H1 H2\"\n"\
" --pulsar            name of pulsar e.g. J0534+2200\n"\
" --par-file          pulsar parameter (.par) file (full path) \n"\
" --input-dir         directory containing the input data files\n"\
" --output-dir        directory for output data files\n"\
" --dob-ul            (double) percentage degree-of-belief for upper limit\n\
                     - if 0 (default no upper limit is produced\n"\
" --output-post       output the full log(posterior)\n"\
" --chunk-min         (int) minimum stationary length of data to be used in\n\
                     the likelihood e.g. 5 mins\n"\
" --chunk-max         (int) maximum stationary length of data to be used in\n\
                     the likelihood e.g. 30 mins\n"\
"\n"\
" Parameter space grid values:-\n"\
" --minh0             (double) minimum of the h0 grid\n"\
" --maxh0             (double) maximum of the h0 grid\n"\
" --h0steps           (int) number of steps in the h0 grid\n"\
" --minphi0           (double) minimum of the phi0 grid\n"\
" --maxphi0           (double) maximum of the phi0 grid\n"\
" --phi0steps         (int) number of steps in the phi0 grid\n"\
" --minpsi            (double) minimum of the psi grid\n"\
" --maxpsi            (double) maximum of the psi grid\n"\
" --psisteps          (int) number of steps in the psi grid\n"\
" --minci             (double) minimum of the cos(iota) grid\n"\
" --maxci             (double) maximum of the cos(iota) grid\n"\
" --cisteps           (int) number of steps in the cos(iota) grid\n"\
" --psi-bins          (int) no. of psi bins in the time-psi lookup table\n"\
" --time-bins         (int) no. of time bins in the time-psi lookup table\n"\
"\n"\
" Prior values:-\n"\
" --h0prior           type of prior on h0 - uniform, jeffreys or gaussian\n"\
" --h0mean            (double) mean of a gaussian prior on h0\n"\
" --h0sig             (double) standard deviation of a gaussian prior on h0\n"\
" --phi0prior         type of prior on phi0 - uniform or gaussian\n"\
" --phi0mean          (double) mean of a gaussian prior on phi0\n"\
" --phi0sig           (double) std. dev. of a gaussian prior on phi0\n"\
" --psiprior          type of prior on psi - uniform or gaussian\n"\
" --psimean           (double) mean of a gaussian prior on psi\n"\
" --psisig            (double) std. dev. of a gaussian prior on psi\n"\
" --ciprior           type of prior on cos(iota) - uniform or gaussian\n"\
" --cimean            (double) mean of a gaussian prior on cos(iota)\n"\
" --cisig             (double) std. dev. of a gaussian prior on cos(iota)\n"\
"\n"\
" MCMC parameters:-\n"\
" --mcmc              1 to perform an MCMC, 0 to not perform an MCMC\n"\
" --interations       (int) the number of iteraction in the MCMC chain\n"\
" --burn-in           (int) the number of burn in iterations\n"\
" --temperature       (double) the temperatue to start of the simulated\n\
                     annealing in the burn in stage\n"\
" --h0-width          (double) width of the h0 proposal distribution (if set\n\
                     to 0 this will be worked out in the code)\n"\
" --psi-width         (double) width of the psi proposal distribution\n"\
" --phi0-width        (double) width of the phi proposal distribution\n"\
" --ci-width          (double) width of the cos(iota) proposal distribution\n"\
" --glitch-times      (char) a string of pulsar glitch times given in MJD\n\
                     e.g. \"45623.872 52839.243 53992.091\" - at each glitch\n\
                     an additional MCMC phase parameter will be used\n"\
" --glitch-cut        (double) the number of seconds of data to ignore\n\
                     before and after a glitch\n"\
"\n"

#define MAXLENGTH 500000

/** define structures */

/* a structure to hold the four instrinsic pulsar parameters */
typedef struct tagIntrinsicPulsarVariables{
  double h0;         /* GW amplitude */
  double phi0;       /* initial GW phase */
  double ci;         /* cos(inclination angle) */
  double psi;        /* polarisation angle */
  
  /* derived variables to allow quicker computation on the grid */
  /* MUST BE DEFINED PRIOR TO THE LIKELIHOOD CALCULATION */
  double Xplus;      /* 0.5*(1+cos^2(iota)) */
  double Xcross;     /* cos(iota) */
  double Xpsinphi_2; /* 0.5*Xplus*sin(phi) */
  double Xcsinphi_2; /* 0.5*Xcross*sin(phi) */
  double Xpcosphi_2; /* 0.5*Xplus*cos(phi) */
  double Xccosphi_2; /* 0.5*Xcross*cos(phi) */
  
  double Aplus;      /* 0.5*h0*(1+cos^2(iota)) */
  double Across;     /* h0*cos(iota) */
  double Apsinphi_2; /* 0.5*Aplus*sin(phi) */
  double Acsinphi_2; /* 0.5*Across*sin(phi) */
  double Apcosphi_2; /* 0.5*Aplus*cos(phi) */
  double Accosphi_2; /* 0.5*Across*cos(phi) */
}IntrinsicPulsarVariables;

typedef struct tagMeshGrid{
  IntrinsicPulsarVariables minVals;
  IntrinsicPulsarVariables maxVals;
  IntrinsicPulsarVariables delta; /* parameter step sizes */
  
  int h0Steps, phiSteps, psiSteps, ciotaSteps;
  
  /* the number of steps in the range for the lookup table */
  int psiRangeSteps;
  int timeRangeSteps;
  
  LALSource pulsar; /* pulsar position values */
  LALDetector detector;
}MeshGrid;

/* a structure to hold the MCMC initial parameters */
typedef struct tagMCMCParams{
  /* MCMC parameters */
  int doMCMC;                      /* flag can be true 1 or false 0 */
  int iterations;                  /* length of MCMC chain */
  int burnIn;                      /* length of burn in period */
  int temperature;                 /* temperature of simulated annealing during
                                      burn in */
  IntrinsicPulsarVariables sigmas; /* standard deviations of the proposal
                                      distributions */
                                      
  char glitchTimes[500];           /* string containing the times of glitches
                                      in MJD */
  int nGlitches;                   /* the number of glitches */
  double glitchCut;                /* the number of seconds of data to ignore
                                      around the time of the glitch */
}MCMCParams;

typedef struct tagPriorVals{
  char *h0Prior, *phiPrior, *psiPrior, *ciotaPrior; /* prior type e.g.
"uniform", "jeffreys" or "gaussian" */

  IntrinsicPulsarVariables vars;

  /* for gaussian prior pass in the mean and standard deviation */
  double meanh0, stdh0;
  double meanphi, stdphi;
  double meanpsi, stdpsi;
  double meanciota, stdciota;
}PriorVals;

/* a structure to hold the command line input values */
typedef struct tagInputParams{
  char detectors[40];
  char pulsar[12];
  char parFile[256];
  
  char inputDir[256];
  char outputDir[256];
  
  /* prior values */
  PriorVals priors;
  
  /* grid values */
  MeshGrid mesh;
  
  /* maximum and minumum length of time over which we assume stationarity of
     the data e.g. min 5 mins and max 30 mins */
  int chunkMin;
  int chunkMax;
  
  /* parameters for an MCMC */
  MCMCParams mcmc;
  
  double dob;                      /* degree of belief for upper limit */
  
  int outputPost; /* flag for whether of not to output the full posterior */
                                      
  int verbose;
}InputParams;


typedef struct tagDataStructure{
  COMPLEX16Vector *data;     /* real and imaginary B_ks */
  REAL8Vector *times;        /* time stamp for each B_k */
  INT4Vector *chunkLengths;  /* length of each data chunk e.g. 30 mins */
  int chunkMax;              /* maximum chunk length to be used e.g. 30 mins */
  int chunkMin;              /* minumum chunk length to be used e.g. 5 mins */

  REAL8Vector *sumData;      /* sum over the data for each chunk (Re(B)^2 +
                               Im(B)^2)*/
  REAL8Vector *sumModel;     /* sum over the model for each chunk (Re(y)^2 +
                               Im(y)^2) not including the h0 scaling */
  REAL8Vector *sumDataModel; /* sum over the model and data for each chunk
                               (Re(B)*Re(y) + Im(B)*Im(y)) */
}DataStructure;

/* a detector response function lookup table structure */
typedef struct tagDetRespLookupTable{
  LALDetAMResponse **lookupTable; /* array containing the lookup table */
  int psiSteps;                   /* number of steps across psi range */
  int timeSteps;                  /* number of steps across a day */
}DetRespLookupTable;

typedef struct tagOutputParams{
  char *outputDir;  /* directory to output files to */
  char *psr;        /* pulsar name */
  char *det;        /* detector e.g. H1, H2, L1, G1, V1, Joint */
  
  char *margParam;  /* parameter over which to not marginalise */
  
  int outPost;      /* flag for whether to output the full posterior 1=yes,
                       0=no */
  double dob;       /* the degree-of-belief for an upper limits calculation -
                       don't perform calc if set to zero */ 
}OutputParams;


/** define functions */

/* function to get the input arguments from the command line */
void get_input_args(InputParams *inputParams, int argc, char *argv[]);

/* function to allocate memory for a likelihood array - it will be index as
   logLike[phi][cosiota][psi][h0] */
double **** allocate_likelihood_memory(MeshGrid mesh);

/* function to create a log likelihood array over the parameter grid */
void create_likelihood_grid(DataStructure data, double ****logLike, 
  MeshGrid mesh);

/* a function to compute the log likelihood of the data given some parameters 
- this function loops over the h0 part of the likelihood array for speed */
void log_likelihood(double *likeArray, DataStructure data,
  IntrinsicPulsarVariables vars, MeshGrid mesh, DetRespLookupTable lookupTable);

/* a function to sum over the data */
void sum_data(DataStructure data);

/* a function to give the log prior given a set of prior ranges and parameter
values */
double log_prior(PriorVals priors, MeshGrid mesh);

/* function to compute the log posterior */
double log_posterior(double ****logLike, PriorVals prior, MeshGrid mesh,
  OutputParams output);

/* marginalise posterior over requested parameter and output the log evidence if
requested */
double marginalise_posterior(double ****logPost, PriorVals prior, 
  MeshGrid mesh, OutputParams output);
  
/* detector response lookup table function  - this function will output a lookup
table of points in time and psi, covering a day from the start time (t0) and
from -pi/4 to pi/4 in psi */
void response_lookup_table(double t0, LALDetAndSource detAndSource,
  DetRespLookupTable lookupTable);

/* function to do the trapezium rule for integration on logs  */
double log_trapezium(double logHeight1, double logHeight2, double width);
  
/* factorial function */
int factorial(int num);

/* function to combine log likelihoods to give a joint likelihood */
void combine_likelihoods(double ****logLike1, double ****logLike2, 
  MeshGrid mesh);
  
#ifdef  __cplusplus
}
#endif

#endif /* _PULSAR_PARAMETER_ESTIMATION_H */
