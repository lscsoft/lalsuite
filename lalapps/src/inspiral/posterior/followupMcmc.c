/*
*  Copyright (C) 2009 Christian Roever
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

/*************************************************/
/*                                               */
/*  Coherent follow-up MCMC code                 */
/*                                               */
/*  Christian Roever                             */
/*                                               */
/* * * * * * * * * * * * * * * * * * * * * * * * */
/*                                               */
/*  Compilation works (for me) using the         */
/*  following command:                           */
/*                                               */
/*  gcc -ansi -O3 -o followupMcmc followupMcmc.c -static -I ${LSCSOFT_LOCATION}/libframe/include -L ${LSCSOFT_LOCATION}/libframe/lib -I ${LSCSOFT_LOCATION}/lal/include -L ${LSCSOFT_LOCATION}/lal/lib -llal -lFrame -lfftw3 -lfftw3f -lgsl -lgslcblas -lm   */
/*                                               */
/* * * * * * * * * * * * * * * * * * * * * * * * */
/*                                               */
/*  Things you may want to play around with      */
/*  are in particular proposals (look out        */
/*  for occurrences of the "setstdev()" and      */
/*  "setcor()" functions, as well as the         */
/*  "proposeInspiralNospin()" function           */
/*  itself).                                     */
/*                                               */
/*************************************************/



#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>   /* FFT library: www.fftw.org                                      */
#include <FrameL.h>  /* Frame library: http://lappweb.in2p3.fr/virgo/FrameL            */
#include <time.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/Units.h>
#include <lal/LALInspiral.h>
#include <gsl/gsl_rng.h>     /* http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html */
#include <gsl/gsl_randist.h> /* http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Distributions.html */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>


/*-- define some physical constants --*/
const double Msun = 1.9889194662e30;   /* solar mass kg */
const double Mpc  = 3.08568025e22;     /* metres in a Mpc (LAL:3.0856775807e22) */
const double G    = 6.67259e-11;       /* 6.674215e-11; */ 
const double c    = 299792458.0;       /* speed of light m/s */ 
const double pi   = 3.141592653589793; 
const double earthRadiusEquator = 6378137.0;           /* (WGS84 value)                    */
const double earthFlattening    = 0.00335281066474748; /* (WGS84 value: 1.0/298.257223563) */
/*const double earthRadiusPole    = 6356752.314; */    /* (WGS84 value)                    */

/*-- set verbosity flag: --*/
const int verbose = 1;

/*-- set downsampling factor (=1 for no downsampling; 2, 4, or 8 otherwise): --*/
const int downsamplingfactor = 8;

/*-- set flag forcing signal into flat part of Tukey window: --*/
const int forceFlatTukey = 1;

/*-- signal "families" available (common parameter sets): --*/
enum signal {
  InspiralNoSpin,    /*  0) 9 parameters  */
  BurstSineGaussian  /*  1) 8 parameters  */
};

/*-- different templates available: --*/
enum template {
  iR2PN,          /*  0) restricted 2.0PN inspiral              */
  i20SP,          /*  1) 2.0PN stationary phase inspiral        */
  i25SP,          /*  2) 2.5PN stationary phase inspiral        */
  i2025,          /*  3) 2.0PN amplitude, 2.5PN phase inspiral  */
  i2535,          /*  4) 2.5PN amplitude, 3.5PN phase inspiral  */
  iLALTT2PN00,    /*  5) LAL Taylor T2 Newtonian                */
  iLALTT2PN10,    /*  6) LAL Taylor T2 1.0 PN                   */
  iLALTT2PN15,    /*  7) LAL Taylor T2 1.5 PN                   */
  iLALTT2PN20,    /*  8) LAL Taylor T2 2.0 PN                   */
  iLALTT3PN00,    /*  9) LAL Taylor T3 Newtonian                */
  iLALTT3PN10,    /* 10) LAL Taylor T3 1.0 PN                   */
  iLALTT3PN15,    /* 11) LAL Taylor T3 1.5 PN                   */
  iLALTT3PN20,    /* 12) LAL Taylor T3 2.0 PN                   */
  iLALIMRPhenomA, /* 13) LAL Phenomenological                   */
  iLALEOBNR,      /* 14) LAL EOBNR                              */
  bSineGaussian   /* 15) sine-gaussian burst                    */
};


/* "interferometer" structure:                                                                  */
typedef struct {
  char   *name;             /* interferometer's name                                            */
  double latitude;          /* (radian)                                                         */
  double longitude;         /* (radian)                                                         */
  double rightArmAngle;     /* right arm direction (counter-clockwise (!) from North, radian)   */
  double positionVector[3]; /* vector pointing from geocentre to ifo, in units of metres (!)    */
  double rightArmVector[3]; /* unit vector pointing along right arm                             */
  double leftArmVector[3];  /* unit vector pointing along left arm                              */
  double normalVector[3];   /* unit normal vector of ifo plane                                  */
} interferometer;


/* "DataFramework" structure:                                                                   */
typedef struct {
  /* model details:                                                                             */
  double       timeCenter;      /* mean of coalecence time's (uniform) prior (GPS seconds)      */
  double       timeBefore;      /* amount of data (seconds) to be read before (prior) t_c       */
  double       timeAfter;       /* amount of data (seconds) to be read after (prior) t_c        */
  double       minF, maxF;      /* Fourier frequency range of interest (Hz)                     */
  double       tukeypar;        /* value of (Tukey-) windowing parameter (0 <= a <= 1)          */
  /* in/output details:                                                                         */
  char         *frameChannel;   /* channel name         (e.g.: "L1:LSC-Strain")                 */
  char         *datacachefile;  /* cache file containing GWF file names                         */
  char         *noisecachefile; /* ditto                                                        */
  char         *datafilenames;  /* string of GWF file names, assembled during initialisation    */
  char         *noisefilenames; /* ditto                                                        */
  double       PsdEstimateStart;/* start time of data stretch to use for spectrum estimation    */
  double       PsdEstimateEnd;  /* end time of data stretch to use for spectrum estimation      */
  int          PsdEstimateN;    /* number of segments averaged over for spectrum estimation     */
  int          simulateData;    /* flag indicating whether to simulate data (or read from file) */
  int          simulatePsd;     /* flag indicating whether to copy PSD or simulate estimation   */
  /* elements referring to data:                                                                */
  double       dataStart;       /* time point of first sample in the data (GPS seconds)         */
  double       dataDeltaT;      /* time resolution of the data (aka `cadence', `dt')            */
  long         dataSize;        /* number of samples in data                                    */
  double       *data;           /* data vector (in time domain)                                 */
  double       *window;         /* (values of) windowing function                               */
  double       winss;           /* sum of squared windowing coefficients                        */
  double complex *dataFT;       /* Fourier transformed (and windowed) data                      */
  long         FTSize;          /* size of `dataFT' vector                                      */
  double       FTDeltaF;        /* frequency resolution of data                                 */
  double       *powspec;        /* logarithmic 1-sided (!) power spectral density values        */
  char         *noiseModel;     /* noise model to be used for simulating noise                  */
  long         minInd, maxInd;  /* range of indices corresponding to "minF" & "maxF"            */
  /* internal variables:                                                                        */
  double       *FTin;           /* (internal) FT input vector (of length `dataSize')            */
  fftw_complex *FTout;          /* (internal) FT output vector (of length `FTSize')             */
  fftw_plan    FTplan;          /* Transform plan (only used internally--access via `FTexec()') */
  double       rawDataRange;    /* raw data to be read to get correct range after downsampling  */
  interferometer *ifo;          /* interferometer details                                       */
} DataFramework;


/* "vector" structure: */
typedef struct {
  int    dimension;
  double *value;
  char   **name;
} vector;


/* "McmcFramework" structure:                                                                     */
typedef struct {
  int    template;           /* index of template to be used                                      */
  int    parameterset;       /* index of (implied) parameter set                                  */
  int    pardim;             /* dimension of parameter space                                      */
  vector fixed;              /* (logical) vector indicating which parameters have fixed values    */
  vector priorparameters;    /* parameters of prior distn. (depends on model used)                */
  vector startvalue;         /* starting value parameter vector                                   */
  double GMST;               /* reference Greenwich mean sidereal time (in radian)                */
  int    celestialCoords;    /* output celestial instead of geographical coordinates              */
  double *covariance;        /* proposal covariance matrix                                        */
  gsl_matrix *covcholesky;   /* Cholesky decomposition of above covariance matrix                 */
  double studentDF;          /* degrees-of-freedom parameter for Student-t proposals (0.0=Normal) */
  double *guessparameters;   /* parameters of "guessing distribution"                             */
  char   *logfilename;       /* name of output log-file                                           */
  long   iterations;         /* number of MCMC iterations to run                                  */
  double secPerIteration;    /* (estimated) seconds per iteration                                 */
} McmcFramework;

gsl_rng *GSLrandom;             /* GSL random number generator */

void vectorInit(vector *vec);
void vectorAdd(vector *vec, char name[], double value);
void vectorSetup(vector *vec, int signal);
void vectorDispose(vector *vec);
void vectorCopy(vector *vec1, vector *vec2);
void vectorSetValue(vector *vec, char name[], double value);
double vectorGetValue(vector *vec, char name[]);
int vectorIsElement(vector *vec, char name[]);
void vectorPrint(vector *vec);

int char2template(char *templatename);

void expandLocationString(char *shortForm, char *longForm);
interferometer *ifoPointer(interferometer *ifolist, int ifoN, char locationString[]);

void parseParameterOptionString(char *input, char **parnames[], double **parvalues, int *n);
void parseCharacterOptionString(char *input, char **strings[], int *n);
int noVectorStrg(char *input);
void printhelpmessage();
int init(DataFramework *DFarg[], McmcFramework *MFarg[],
         int *coherentN,
         int argc, char *argv[],
         interferometer *ifolist, int ifoN);
void clearDF(DataFramework *DF);
int readData(DataFramework *DF);
char *cache2string(char *cachefilename, double *startTime, double *endTime, char *locationstrg);
void generateNoise(double deltat, int N, 
                   double* noise, char *model);
void simulateData(DataFramework *DF);
int estimatePSD(DataFramework *DF);
void simulatePsdEstimation(DataFramework *DF);
double *downsample(double data[], int *datalength, int factor);
double tukeywindow(int j, int N, double r);
void FTexec(DataFramework *DF, double *input, fftw_complex *output);
double spectrumIniLigo(double f);
double spectrumAdvLigo(double f);
double spectrumVirgo(double f);
double spectrumGeo(double f);
void normalise(double vec[3]);
void rotate(double x[3], double angle, double axis[3]);
int righthanded(double x[3], double y[3], double z[3]);
void orthoproject(double x[3], double vec1[3], double vec2[3]);
double angle(double x[3], double y[3]);
void coord2vec(double lati, double longi, double x[3]);
void vec2coord(double x[3], double *lati, double *longi);
void ifoInitComputations(interferometer *ifolist, int ifoN);
void ifoInit(interferometer *ifolist[], int *ifoN);
void clearIfo(interferometer *ifolist, int *ifoN);
void localParameters(vector *parameter, interferometer *ifo,
                     double *timeshift, double *polarisation, double *altitude, double *azimuth);

void antennaepattern(double altitude, double azimuth, double polarisation,
                     double *Fplus, double *Fcross);
double mc2mass1(double mc, double eta);
double mc2mass2(double mc, double eta);
double mc2mt(double mc, double eta);
double logJacobianMcEta(double mc, double eta);

void signaltemplate(DataFramework *DF, int waveform, vector *parameter, double complex *output);
void template2025(DataFramework *DF, vector *parameter, double Fplus, double Fcross,
                  double complex *output);
void template2535(DataFramework *DF, vector *parameter, double Fplus, double Fcross,
                  double complex *output);
void templateStatPhase(DataFramework *DF, vector *parameter, double Fplus, double Fcross,
                       double complex *output, double PNOrder);
void templateR2PN(DataFramework *DF, vector *parameter, double Fplus, double Fcross,
                  double complex *output);
void templateLAL(DataFramework *DF, vector *parameter, double Fplus, double Fcross,
                 double complex *output, int approximant, int order);
void templateSineGaussianBurst(DataFramework *DF, vector *parameter, double Fplus, double Fcross,
                               double complex *output);
void inject(DataFramework *DF, int coherentN, int waveform, vector *parameter);
void dumptemplates(DataFramework *DF, vector *parameter, char *filenameF, char *filenameT);

double loglikelihood(DataFramework *DF, int coherentN, int waveform, vector *parameter);
double signaltonoiseratio(DataFramework *DF, int coherentN, int waveform, vector *parameter,
                          double indivSNRs[]);

void clearMF(McmcFramework *MF);

void setcov(McmcFramework *MF, int parameter1, int parameter2, double covariance);
void setvar(McmcFramework *MF, int parameter, double variance);
void setstdev(McmcFramework *MF, int parameter, double standarddeviation);
double getcov(McmcFramework *MF, int parameter1, int parameter2);
double getvar(McmcFramework *MF, int parameter);
double getstdev(McmcFramework *MF, int parameter);
void setcor(McmcFramework *MF, int parameter1, int parameter2, double correlation);
void RandMVNorm(McmcFramework *MF, double *result);

double logprior(McmcFramework *MF, vector *parameter);
double logpriorInspiralNospin(McmcFramework *MF, vector *parameter);
double logpriorSineGaussianBurst(McmcFramework *MF, vector *parameter);

double logguessdensity(McmcFramework *MF, vector *parameter);
double logguessInspiralNospin(McmcFramework *MF, vector *parameter);

void priordraw(McmcFramework *MF, vector *parameter);
void priordrawInspiralNospin(McmcFramework *MF, vector *parameter);
void priordrawBurstSineGaussian(McmcFramework *MF, vector *parameter);

void guess(McmcFramework *MF, vector *parameter);
void guessInspiralNospin(McmcFramework *MF, vector *parameter);

void propose(McmcFramework *MF, DataFramework *DF, int coherentN, 
             vector *parameter, double *logMHcoef);
void proposeInspiralNospin(McmcFramework *MF, DataFramework *DF, int coherentN, 
                           vector *parameter, double *logMHcoef);
void proposeBurstSineGaussian(McmcFramework *MF, DataFramework *DF, int coherentN, 
                              vector *parameter, double *logMHcoef);

void importanceresample(DataFramework *DF, int coherentN, McmcFramework *MF,
                        vector *parameter, 
                        long samplesize, long subsamplesize);
void logtofile(McmcFramework *MF, vector *parameter, 
               long iteration, long accepted, 
               double logprior, double loglikelihood, double logposterior);
void printtime();
void savePSD(DataFramework *DF, char *filename);
void metropolishastings(McmcFramework *MF, DataFramework *DF, int coherentN);
void printDF(DataFramework *DF);
int numberIfoSites(DataFramework *DF, int coherentN);
double GMST(double GPSsec);
double rightAscension(double longi, double gmst);
double longitude(double rightascension, double gmst);

void vectorInit(vector *vec)
/* Very basic initialisation of a vector. */
/* NEEDS TO BE RUN FOR ANY NEW VECTOR !!  */
{
  vec->dimension = 0;
  vec->value     = NULL;
  vec->name      = NULL;
}


void vectorAdd(vector *vec, char name[], double value)
/* Takes a given vector and adds a new dimension (name & value). */
{
  int newdim=0;
  double *newvalue=NULL;
  char **newname=NULL;
  int i,j,k;
  /* allocate new pointers etc: */
  newdim   = vec->dimension + 1;
  newvalue = (double*) malloc(sizeof(double)*newdim);
  newname  = (char**) malloc(sizeof(char*)*newdim);
  if (name[0]=='\0') printf(" : WARNING: empty vector name in 'vectorAdd(\"\",%.3e)'!\n",value);
  if (vec->dimension > 0){
    /* first check for duplicate name: */
    j=0;
    for (i=0; i<vec->dimension; ++i)
      j += (strcmp(vec->name[i], name)==0);
    if (j!=0) printf(" : WARNING: duplicate vector name in 'vectorAdd(\"%s\",%.3e)'!\n",
                     name, value);
    /* copy previous stuff over: */
    for (i=0; i<vec->dimension; ++i) {
      newvalue[i] = vec->value[i];
      j = 0;
      while (vec->name[i][j] != '\0') ++j;
      /* ('j' is number of non-zero characters) */
      newname[i] = (char*) malloc(sizeof(char) * (j+1));
      for (k=0; k<j; ++k)
        newname[i][k] = vec->name[i][k];
      newname[i][j] = '\0';
    }
    /* dispose old stuff: */
    free(vec->value);
    for (i=0; i<vec->dimension; ++i)
      free(vec->name[i]);
    free(vec->name);
  }
  /* add new value & name: */
  newvalue[newdim-1] = value;
  j=0;
  while (name[j] != '\0') ++j;
  newname[newdim-1] = (char*) malloc(sizeof(char)*(j+1));
  for (k=0; k<j; ++k)
    newname[newdim-1][k] = name[k];
  newname[newdim-1][j] = '\0';
  /* link new values & names to old vector: */
  vec->dimension = newdim;
  vec->value     = newvalue;
  vec->name      = newname;
}


void vectorSetup(vector *vec, int signal)
/* Set up parameter vector according to given signal type. */
/* !! Vector needs to be pre-initialised !!                */
{
  if (vec->dimension > 0){
    printf(" : ERROR: attempt to re-setup an existing parameter vector in 'vectorSetup()'!\n");
  }
  else{
    /* parameters used for ANY signal: */
    vectorAdd(vec, "latitude", 0.0);      /* w.r.t. earth frame (geographic lat.)  */
    vectorAdd(vec, "longitude", 0.0);     /* w.r.t. earth frame (geographic long.) */
    vectorAdd(vec, "polarisation", 0.0);  /* w.r.t. earth frame                    */
    vectorAdd(vec, "time", 0.0);          /* w.r.t. earth frame (at geocentre)     */
    /* signal-specific parameters: */
    if (signal == InspiralNoSpin){
      vectorAdd(vec, "phase", 0.0);
      vectorAdd(vec, "logdistance", 0.0);
      vectorAdd(vec, "chirpmass", 0.0);
      vectorAdd(vec, "massratio", 0.0);
      vectorAdd(vec, "inclination", 0.0);
    }
    else if (signal == BurstSineGaussian){
      vectorAdd(vec, "phase", 0.0);
      vectorAdd(vec, "logamplitude", 0.0);
      vectorAdd(vec, "logsigma", 0.0);
      vectorAdd(vec, "frequency", 0.0);
    }
    else 
      printf(" : ERROR: unknown signal type specification in 'vectorSetup()'!\n");
  }
}


void vectorDispose(vector *vec)
/* Wipes out a vector's contents and frees memory.         */
/* A 'disposed' vector does not need to be re-initialised. */
{
  int i;
  if (vec->dimension > 0) {
    free(vec->value);
    for (i=0; i<vec->dimension; ++i)
      free(vec->name[i]);
    free(vec->name);
  }
  vec->dimension = 0;
  vec->value     = NULL;
  vec->name      = NULL;
}


void vectorCopy(vector *vec1, vector *vec2)
/* Copy vec1 over to (pre-initialised!!) vec2 */
{
  int i,j,k;
  if (vec2->dimension != vec1->dimension) {
    vectorDispose(vec2);
    vec2->dimension = vec1->dimension;
    if (vec1->dimension > 0) {
      vec2->value = (double*) malloc(sizeof(double)*vec1->dimension);
      vec2->name  = (char**) malloc(sizeof(char*)*vec1->dimension);
    }
  }
  if (vec1->dimension > 0) {
    for (i=0; i<vec1->dimension; ++i){
      vec2->value[i] = vec1->value[i];
      j=0;
      while (vec1->name[i][j] != '\0') ++j;
      vec2->name[i] = (char*) malloc(sizeof(char)*(j+1));
      for (k=0; k<j; ++k)
        vec2->name[i][k] = vec1->name[i][k];
      vec2->name[i][j] = '\0';
    }
  }
  else {
    vec2->value = NULL;
    vec2->name  = NULL;
  }
}


void vectorSetValue(vector *vec, char name[], double value)
{
  int i, notfound=1;
  if (vec->dimension > 0) {
    i=0;
    while ((i<vec->dimension) && (notfound=(strcmp(vec->name[i], name)!=0))) ++i;
    if (notfound)
      printf(" : ERROR: attempt to set unknown vector element in 'vectorGetValue(...,\"%s\")'!\n", 
             name);
    else
      vec->value[i] = value;
  }
  else
    printf(" : ERROR: attempt to set element of empty vector in 'vectorGetValue(...,\"%s\")'!\n", 
           name);
}


double vectorGetValue(vector *vec, char name[])
{
  int i, notfound=1;
  double result=0.0;
  if (vec->dimension > 0) {
    i=0;
    while ((i<vec->dimension) && (notfound=(strcmp(vec->name[i], name)!=0))) ++i;
    if (notfound)
      printf(" : ERROR: requesting unknown vector element in 'vectorGetValue(...,\"%s\")'!\n", 
             name);
    else
      result = vec->value[i];
  }
  else
    printf(" : ERROR: requesting element from empty vector in 'vectorGetValue(...,\"%s\")'!\n", 
           name);
  return result;
}


int vectorIsElement(vector *vec, char name[])
/* returns 1 if 'vec' has an element named 'name', and zero otherwise. */
{
  int i, notfound=1;
  int result = 0;
  if (vec->dimension > 0) {
    i=0;
    while ((i<vec->dimension) && (notfound=(strcmp(vec->name[i], name)!=0))) ++i;
    if (!notfound)
      result = 1;
  }
  return result;
}


void vectorPrint(vector *vec)
{
  int i;
  printf(" : <vector(%d):", vec->dimension);
  if (vec->dimension == 0)
    printf(" EMPTY");
  else {
    for (i=0; i<vec->dimension; ++i) {
      printf(" '%s'=%.3e", vec->name[i], vec->value[i]);
    }
  }
  printf(">\n");
}


int char2template(char *templatename)
/* convert character strings (as supplied in command line arguments) to integers. */
{
  int result = -1;
  if (strcmp(templatename, "R2PN")==0)  
    result = iR2PN;
  else if (strcmp(templatename, "20SP")==0)  
    result = i20SP;
  else if (strcmp(templatename, "25SP")==0)  
    result = i25SP;
  else if (strcmp(templatename, "2025")==0)  
    result = i2025;
  else if (strcmp(templatename, "2535")==0)  
    result = i2025;

  else if (strcmp(templatename, "LALTaylorT2PN00")==0)  
    result = iLALTT2PN00;
  else if (strcmp(templatename, "LALTaylorT2PN10")==0)  
    result = iLALTT2PN10;
  else if (strcmp(templatename, "LALTaylorT2PN15")==0)  
    result = iLALTT2PN15;
  else if (strcmp(templatename, "LALTaylorT2PN20")==0)  
    result = iLALTT2PN20;
  else if (strcmp(templatename, "LALTaylorT3PN00")==0)  
    result = iLALTT3PN00;
  else if (strcmp(templatename, "LALTaylorT3PN10")==0)  
    result = iLALTT3PN10;
  else if (strcmp(templatename, "LALTaylorT3PN15")==0)  
    result = iLALTT3PN15;
  else if (strcmp(templatename, "LALTaylorT3PN20")==0)  
    result = iLALTT3PN20;
  else if (strcmp(templatename, "LALIMRPhenomA")==0) {
    result = iLALIMRPhenomA;
    printf(" :\n : WARNING: Phenomenological templates do not yet work !!\n :\n");
  }
  else if (strcmp(templatename, "LALEOBNR")==0)  
    result = iLALEOBNR;

  else if (strcmp(templatename, "SineGaussian")==0)  
    result = bSineGaussian;
  return result;
}


void expandLocationString(char *shortForm, char *longForm)
/* expands a "short" detector location identification string                    */
/* (as used in cache files and for the "--network" command line argument)       */
/* to the "long" form that can be found in the "interferometer.name" slot.      */
/* Short forms: "H", "L", "V", "G"                                              */
/* long forms : "LIGO-Hanford", "LIGO-Livingston", "VIRGO-Pisa", "GEO-Hannover" */
{
  if (strcmp(shortForm, "H")==0)
    strcpy(longForm, "LIGO-Hanford");
  else if (strcmp(shortForm, "L")==0)
    strcpy(longForm, "LIGO-Livingston");
  else if (strcmp(shortForm, "V")==0)
    strcpy(longForm, "VIRGO-Pisa");
  else if (strcmp(shortForm, "G")==0)
    strcpy(longForm, "GEO-Hannover");
  else {
    strcpy(longForm, "??");
    printf(" : ERROR in 'expandLocationString()': unknown abbreviation \"%s\".\n", shortForm);
  }
}


interferometer *ifoPointer(interferometer *ifolist, int ifoN, char locationString[])
/* for a given list of interferometers "ifolist" of length "ifoN",     */
/* the list is searched for an interferometer named "locationString",  */
/* and a pointer to the corresponding list entry is returned.          */
{
  interferometer *result = NULL;
  int found = 0;
  int i = 0;
  while ((!found) && (i<ifoN)) {
    found = (strcmp(locationString, ifolist[i].name)==0);
    if (!found) ++i;
  }
  result = &ifolist[i];
  if (!found) {
    printf(" : ERROR in 'ifoPointer()': unknown location \"%s\".\n", locationString);
    result = NULL;
  }
  return result;
}


void parseParameterOptionString(char *input, char **parnames[], double **parvalues, int *n)
/* parses a character string (passed as one of the options) and decomposes */
/* it into parameter names and -values. Input is of the form               */
/*   input    :  "[mc=3.1415,eta=0.24,tc=700009616.123]"                   */
/* and the resulting output is                                             */
/*   parnames :  {"mc", "eta", "tc"}                                       */
/*   parvalues:  {3.1415, 0.24, 700009616.123}                             */
/* length of parameter names is by now limited to 31 characters.           */
{
  int i,j,k,l,m;
  char strg[32];

  /* perform a very basic well-formedness-check and count number of parameters: */
  i = j = 0;
  *n = 0;
  while (input[i] != '\0') {
    if ((j==0) & (input[i] == '[')) j = 1;
    if ((j==1) & (input[i] == '=')) j=2;
    if ((j==1) & (input[i] == ',')) {j=1; ++(*n);}
    if ((j==1) & (input[i] == ']')) {j=3; ++(*n);}
    if ((j==2) & (input[i] == ',')) {j=1; ++(*n);}
    if ((j==2) & (input[i] == ']')) {j=3; ++(*n);}
    ++i;
  }
  if (j!=3) printf(" : ERROR: argument vector '%s' not well-formed!\n", input);

  /* allocate memory for results: */
  *parvalues = (double*) malloc(sizeof(double) * (*n));
  *parnames  = (char**)  malloc(sizeof(char*) * (*n));
  for (i=0; i<(*n); ++i) (*parnames)[i] = (char*) malloc(sizeof(char)*32);

  /* copy elements: */
  i = j = l = 0;
  while ((input[i] != '\0')) {
    if ((j==0) & (input[i] == '[')) {j=1; k=i+1;}
    if ((j==1) & (input[i] == '=')) { /* copy preceding variable name */
      j=2;
      /* copy elements  k to i-1 ... */
      for (m=0; m<(i-k); ++m)
        (*parnames)[l][m] = input[k+m];
      (*parnames)[l][m] = '\0';
      k=i+1;
    }
    if ((j==1) & ((input[i] == ',') | (input[i] == ']'))) { /* no preceding variable name, copy value */
      j=1;
      (*parnames)[l][0] = '\0';
      for (m=0; m<(i-k); ++m)
        strg[m] = input[k+m];
      strg[m] = '\0';
      (*parvalues)[l] = atof(strg);
      k=i+1;
      ++l;
    }
    if ((j==2) & ((input[i] == ',')|(input[i] == ']'))) { /* copy preceding value */
      j=1; 
      /* copy elements  k to i-1 ... */
      for (m=0; m<(i-k); ++m)
        strg[m] = input[k+m];
      strg[m] = '\0';
      (*parvalues)[l] = atof(strg);
      k=i+1;
      ++l;
    }
    ++i;
  }
}


void parseCharacterOptionString(char *input, char **strings[], int *n)
/* parses a character string (passed as one of the options) and decomposes   */
/* it into individual parameter character strings. Input is of the form      */
/*   input   :  "[one,two,three]"                                            */
/* and the resulting output is                                               */
/*   strings :  {"one", "two", "three"}                                      */
/* length of parameter names is by now limited to 512 characters.            */
/* (should 'theoretically' (untested) be able to digest white space as well. */
/* Irrelevant for command line options, though.) */
{
  int i,j,k,l;
  /* perform a very basic well-formedness-check and count number of parameters: */
  i=0; j=0;
  *n = 0;
  while (input[i] != '\0') {
    if ((j==0) & (input[i]=='[')) j=1;
    if ((j==1) & (input[i]==',')) ++*n;
    if ((j==1) & (input[i]==']')) {++*n; j=2;}
    ++i;
  }
  if (j!=2) printf(" : ERROR: argument vector '%s' not well-formed!\n", input);

  /* now allocate memory for results: */
  *strings  = (char**)  malloc(sizeof(char*) * (*n));
  for (i=0; i<(*n); ++i) (*strings)[i] = (char*) malloc(sizeof(char)*512);

  i=0; j=0; 
  k=0; /* string counter    */
  l=0; /* character counter */
  while ((input[i] != '\0') & (j<3)) {
    /* state transitions: */
    if ((j==0) & ((input[i]!='[') & (input[i]!=' '))) j=1;
    if (((j==1)|(j==2)) & (input[i]==',')) {(*strings)[k][l]='\0'; j=2; ++k; l=0;}
    if ((j==1) & (input[i]==' ')) j=2;
    if ((j==1) & (input[i]==']')) {(*strings)[k][l]='\0'; j=3;}
    if ((j==2) & (input[i]==']')) {(*strings)[k][l]='\0'; j=3;}
    if ((j==2) & ((input[i]!=']') & (input[i]!=',') & (input[i]!=' '))) j=1;
    /* actual copying: */
    if (j==1) {
      if (l>=511) {
        printf(" : WARNING: character argument too long!\n");
        printf(" : \"%s\"\n",(*strings)[k]);
      }
      else {
        (*strings)[k][l] = input[i];
        ++l;
      }
    }
    ++i;
  } 
}


int noVectorStrg(char *input)
/* this function does a (basic) check whether the input string "input"    */
/* looks like a vector of the form as digested by the above two functions */
/* "parseParameterOptionString()" or "parseCharacterOptionString()".      */
/* (basically checks for opening and closing brackets.)                   */
{
  int i, j;
  int result;
  i=0; j=0;
  while (input[i] != '\0') {
    if ((j==0) & (input[i]=='[')) j=1;
    if ((j==1) & (input[i]==']')) j=2;
    ++i;
  }
  if (j!=2) result = 1;
  else result = 0;
  return result;
}


void printhelpmessage()
{
  printf(" | \n");
  printf(" | Coherent follow-up MCMC code.\n");
  printf(" | \n");
  printf(" | Options:              (default:)         \n");
  printf(" |   --help                                 print this message        \n");
  printf(" |   --template                             specify signal template   \n");
  printf(" |   --logfilename                          output (text) file name   \n");
  printf(" |   --iterations          (1000000)        number of MCMC iterations \n");
  printf(" |   --randomseed                           1-2 integers for initialisation\n");
  printf(" |   --tcenter                              'center' of data to be analysed\n");
  printf(" |   --tbefore             (30)             margin before...               \n");
  printf(" |   --tafter              (1)              ...and after above 'center'    \n");
  printf(" |   --tukey               (0.1)            parameter for Tukey window     \n");
  printf(" |   --cachefile                            cache file for GWF files     \n");
  printf(" |   --network                              interferometer sites\n");
  printf(" |   --filechannel                          channel name for gwf files\n");
  printf(" |   --specdens            ('initialLigo')  spectral density to be used\n");
  printf(" |   --fixspecdens                          fix PSD (instead of estimating)\n");
  printf(" |   --psdestimatestart                     starting time for PSD estimation\n");
  printf(" |   --psdestimateend                       end time for PSD estimation\n");
  printf(" |   --psdestimaten        (100)            number of segments to average over\n");
  printf(" |   --freqlower           (40)             min. frequency considered\n");
  printf(" |   --frequpper           (1800)           max. frequency considered\n");
  printf(" |   --fixed                                vector of fixed parameters\n");
  printf(" |   --start                                vector of starting values\n");
  printf(" |   --guess                                vector of parameter guesses\n");
  printf(" |   --inject                               vector of injection parameters\n");
  printf(" |   --injecttemplate      (--template)     template to be used for injection\n");
  printf(" |   --importanceresample                   number of startvalue draws\n");
  printf(" |   --priorparameters                      vector of prior parameters\n");
  printf(" | \n");
  printf(" | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  printf(" | Example 1 -- Simulate data, inject a signal and recover it:\n");
  printf(" | \n");
  printf(" | ./followupMcmc --network [H,L,V] --randomseed [123,456] --tcenter 100 --logf\n");
  printf(" | ilename /home/user/data/example01.txt --template 25SP --specdens [initialLig\n");
  printf(" | o,initialLigo,Virgo] --freqlower [40,40,30] --tbefore [20.0,20.0,30.0] --inj\n");
  printf(" | ect [chirpmass=2.0,massratio=0.24,inclination=1.0,time=100.0,logdistance=3.2\n");
  printf(" | ,latitude=0.5,longitude=-1.0,polarisation=1.0,phase=3.0] --guess [chirpmass=\n");
  printf(" | 2.0,massratio=0.24,time=100,distance=25.0] --importanceresample 10000\n");
  printf(" | \n");
  printf(" | Example 2 -- Run MCMC on data read from files:\n");
  printf(" | \n");
  printf(" | ./followupMcmc --cachefile [/home/user/data/H1.cache,/home/user/data/H2.cach\n");
  printf(" | e,/home/user/data/L1.cache] --template 25SP --tcenter 873739911.131 --tbefor\n");
  printf(" | e 20.0 --filechannel [H1:LSC-STRAIN,H2:LSC-STRAIN,L1:LSC-STRAIN] --psdestima\n");
  printf(" | testart 873737200 --psdestimateend 873739500 --importanceresample 100000 --r\n");
  printf(" | andomseed 123 --logfilename /home/user/data/example02.txt --iterations 10000\n");
  printf(" | 000 --priorparameters [0.75,7,873739911.081,873739911.181,40,80] --guess [2.\n");
  printf(" | 0,0.13,873739911.131,40]\n");
  printf(" | \n");
  printf(" | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  printf(" | There are currently 13 different waveform templates implemented:\n");
  printf(" |  - 'InspiralNoSpin' templates (9 parameters):\n");
  printf(" |    - 2.0 PN stationary phase approximation          ('20SP')\n");
  printf(" |    - 2.5 PN stationary phase approximation          ('25SP')\n");
  printf(" |    - 2.0 PN amplitude / 2.5 PN phase approximation  ('2025')\n");
  printf(" |    - 2.5 PN amplitude / 3.5 PN phase approximation  ('2535')\n");
  printf(" |    - LAL Taylor T2, Newtonian                       ('LALTaylorT2PN00')\n");
  printf(" |    - LAL Taylor T2, 1.0PN                           ('LALTaylorT2PN10')\n");
  printf(" |    - LAL Taylor T2, 1.5PN                           ('LALTaylorT2PN15')\n");
  printf(" |    - LAL Taylor T2, 2.0PN                           ('LALTaylorT2PN20')\n");
  printf(" |    - LAL Taylor T3, Newtonian                       ('LALTaylorT3PN00')\n");
  printf(" |    - LAL Taylor T3, 1.0PN                           ('LALTaylorT3PN10')\n");
  printf(" |    - LAL Taylor T3, 1.5PN                           ('LALTaylorT3PN15')\n");
  printf(" |    - LAL Taylor T3, 2.0PN                           ('LALTaylorT3PN20')\n");
  printf(" |  - 'BurstSineGaussian' templates (8 parameters):\n");
  printf(" |    - Sine-Gaussian burst template                   ('SineGaussian')\n");
  printf(" | The parameters common to all templates are:\n");
  printf(" |  - 'latitude'     (radian)\n");
  printf(" |  - 'longitude'    (radian)\n");
  printf(" |  - 'polarisation' (radian)\n");
  printf(" |  - 'time'         (GPS seconds)\n");
  printf(" |  - 'phase'        (radian)\n");
  printf(" | The parameters specific to 'InspiralNoSpin' templates are:\n");
  printf(" |  - 'logdistance'  (log(Mpc))\n");
  printf(" |  - 'chirpmass'    (sunmass)\n");
  printf(" |  - 'massratio'    (-)\n");
  printf(" |  - 'inclination'  (radian)\n");
  printf(" | The parameters specific to 'BurstSineGaussian' templates are:\n");
  printf(" |  - 'logamplitude' (log(-))\n");
  printf(" |  - 'logsigma'     (log(seconds))\n");
  printf(" |  - 'frequency'    (Hz)\n");
  printf(" | The templates from the two 'InspiralNoSpin' and 'BurstSineGaussian' classes\n");
  printf(" | by now each share common priors and proposal distributions.\n");
  printf(" | \n");
  printf(" | Some command line arguments are character strings, like the specifications\n");
  printf(" | of templates or noise spectral densities. The waveform identifiers are shown\n");
  printf(" | in the enumeration above. The allowed values for power spectral densities of\n");
  printf(" | noise to be generated are by now 'initialLigo' and 'advancedLigo'.\n");
  printf(" | Some command line arguments are vectors of numbers, see e.g. the example\n");
  printf(" | above. These are in particular the '--fixed', '--start', '--inject',\n");
  printf(" | '--randomseed', '--frequency' and '--priorparameters' arguments. The former\n");
  printf(" | three need to be labelled (like the '--inject' argument in the above\n");
  printf(" | example), the others need to be in correct order.\n");
  printf(" | \n");
  printf(" | The parameters of the 'InspiralNoSpin' prior are (in this order):\n");
  printf(" |  - lower mass bound, upper mass bound\n");
  printf(" |  - lower and upper coalescence time bounds\n");
  printf(" |  - distances (in Mpc) at which a 2+2-Ms inspiral has 90%%\n");
  printf(" |    and 10%% detection probability, respectively.\n");
  printf(" | The parameters of the 'BurstSineGaussian' prior are (in this order):\n");
  printf(" |  - lower frequency bound, upper frequency bound\n");
  printf(" |  - lower and upper time ('mu') bounds\n");
  printf(" |  - expected amplitude\n");
  printf(" |  - expected sigma\n");
  printf(" | \n");
  printf(" | For more details, please see the manual.\n");
  printf(" | \n");
}


int init(DataFramework *DFarg[], McmcFramework *MFarg[],
         int *coherentN,
         int argc, char *argv[],
         interferometer *ifolist, int ifoN)
/* initialisation of "DataFramework" and "McmcFramework" structures, */
/* according to supplied command-line options.                       */
{
  DataFramework *DF = NULL;
  McmcFramework *MF = NULL;
  long i, j, importancedraws;
  unsigned long int seed1=12345, seed2=67890; /* random seeds */
  char CLtemplate[512], CLlogfilename[512], CLiterations[512], CLrandomseed[512], 
       CLtimeCenter[512], CLtimeBefore[512], CLtimeAfter[512], CLtukeypar[512], 
       CLsimulateData[512], CLframeChannel[512],
       CLPsdEstimateStart[512], CLPsdEstimateEnd[512], CLPsdEstimateN[512], 
       CLcachefile[512], CLnetwork[512], CLnoiseModel[512], CLsimulatePsd[512], 
       CLfrequencyLower[512], CLfrequencyUpper[512], 
       CLfixedpar[512], CLstartpar[512], CLguesspar[512], CLinjectpar[512], 
       CLinjecttemplate[512], CLimportancedraws[512], CLpriorpar[512];
  char **argNames;
  double *argValues=NULL;
  int argN;

  int tmpInt;
  double tmpDbl;
  char tmpChar1[64], tmpChar2[64];

  int injecttemplate;
  vector injectpar;
  int InitialisationOK = 1;
  char ifoName[32];
  double startGPS, endGPS;
  static struct option long_options[] = {
    {"template",           required_argument, 0, 't'},
    {"logfilename",        required_argument, 0, 'l'},
    {"iterations",         required_argument, 0, 'i'},
    {"randomseed",         required_argument, 0, 'r'},
    {"tcenter",            required_argument, 0, 'c'},
    {"tbefore",            required_argument, 0, 'b'},
    {"tafter",             required_argument, 0, 'a'},
    {"tukey",              required_argument, 0, 'T'},
    {"filechannel",        required_argument, 0, 'C'},
    {"psdestimatestart",   required_argument, 0, 'R'},
    {"psdestimateend",     required_argument, 0, 'e'},
    {"psdestimaten",       required_argument, 0, 'n'},
    {"cachefile",          required_argument, 0, 's'},
    {"network",            required_argument, 0, 'w'},
    {"specdens",           required_argument, 0, 'd'},
    {"fixspecdens",        no_argument,       0, 'D'},
    {"freqlower",          required_argument, 0, 'f'},
    {"frequpper",          required_argument, 0, 'F'},
    {"fixed",              required_argument, 0, 'x'},
    {"start",              required_argument, 0, 'v'},
    {"guess",              required_argument, 0, 'g'},
    {"inject",             required_argument, 0, 'V'},
    {"injecttemplate",     required_argument, 0, 'I'},
    {"importanceresample", required_argument, 0, 'm'},
    {"priorparameters",    required_argument, 0, 'O'},
    {"help",               no_argument,       0, 'h'},
    {0,                    0,                 0,   0}
  };
  int optparse, optionIndex;

  sprintf(ifoName, "LIGO-Livingston\0");

  /* initialise strings to take command line arguments: */
  CLtemplate[0]=0; CLlogfilename[0]=0; CLiterations[0]=0; CLrandomseed[0]=0; 
  CLtimeCenter[0]=0; CLtimeBefore[0]=0; CLtimeAfter[0]=0; CLtukeypar[0]=0; 
  CLsimulateData[0]=0; CLframeChannel[0]=0;
  CLPsdEstimateStart[0]=0; CLPsdEstimateEnd[0]=0; CLPsdEstimateN[0]=0;
  CLcachefile[0]=0; CLnetwork[0]=0; CLnoiseModel[0]=0; CLsimulatePsd[0]=0; 
  CLfrequencyLower[0]=0; CLfrequencyUpper[0]=0; 
  CLfixedpar[0]=0; CLstartpar[0]=0; CLguesspar[0]=0; CLinjectpar[0]=0; 
  CLinjecttemplate[0]=0; CLimportancedraws[0]=0; CLpriorpar[0]=0;

  /* set default values: */
  strcpy(CLtimeBefore,      "30.0");
  strcpy(CLtimeAfter,       "1.0");
  strcpy(CLtukeypar,        "0.10");
  strcpy(CLfrequencyLower,  "40");
  strcpy(CLfrequencyUpper,  "0");
  strcpy(CLsimulateData,    "1");
  strcpy(CLPsdEstimateN,    "100");
  strcpy(CLsimulatePsd,     "1");
  strcpy(CLimportancedraws, "0");
  strcpy(CLiterations,      "1000000");
  strcpy(CLnoiseModel,      "initialLigo");
  injecttemplate = -1;

  /* read in ALL command line arguments before proceeding further   */
  /* (in particular the number of interferometers needs to be known */
  /* before memory can be properly allocated).                      */

  /* loop over command line arguments: */
  while (1) {
    optparse = getopt_long(argc, argv,
                           "t::l::i::r::c::b::a::T::C::R::e::n::s::w::d::D::f::F::x::v::g::V::I::m::O::h::",
                           long_options, &optionIndex);
    if (optparse == -1) break;
    switch (optparse) {
      case 't': {strcpy(CLtemplate, optarg); break;}            /* --template                    */
      case 'l': {strcpy(CLlogfilename, optarg); break;}         /* --logfilename                 */
      case 'i': {strcpy(CLiterations, optarg); break;}          /* --iterations                  */
      case 'r': {strcpy(CLrandomseed, optarg); break;}          /* --randomseed                  */
      case 'c': {strcpy(CLtimeCenter, optarg); break;}          /* --tcenter                     */
      case 'b': {strcpy(CLtimeBefore, optarg); break;}          /* --tbefore                     */
      case 'a': {strcpy(CLtimeAfter, optarg); break;}           /* --tafter                      */
      case 'T': {strcpy(CLtukeypar, optarg); break;}            /* --tukey                       */
      case 'C': {strcpy(CLframeChannel, optarg); break;}        /* --filechannel                 */
      case 'R': {strcpy(CLPsdEstimateStart, optarg); break;}    /* --psdestimatestart            */
      case 'e': {strcpy(CLPsdEstimateEnd, optarg); break;}      /* --psdestimateend              */
      case 'n': {strcpy(CLPsdEstimateN, optarg); break;}        /* --psdestimaten                */
      case 's': {strcpy(CLsimulateData, "0"); strcpy(CLcachefile,optarg); break;} /* --cachefile */
      case 'w': {strcpy(CLsimulateData, "1"); strcpy(CLnetwork,optarg); break;}   /* --network   */
      case 'd': {strcpy(CLnoiseModel, optarg); break;}          /* --specdens                    */
      case 'D': {strcpy(CLsimulatePsd, "0"); break;}            /* --fixspecdens                 */
      case 'f': {strcpy(CLfrequencyLower, optarg); break;}      /* --freqlower                   */
      case 'F': {strcpy(CLfrequencyUpper, optarg); break;}      /* --frequpper                   */
      case 'x': {strcpy(CLfixedpar, optarg); break;}            /* --fixed                       */
      case 'v': {strcpy(CLstartpar, optarg); break;}            /* --start                       */
      case 'g': {strcpy(CLguesspar, optarg); break;}            /* --guess                       */
      case 'V': {strcpy(CLinjectpar, optarg); break;}           /* --inject                      */
      case 'I': {strcpy(CLinjecttemplate, optarg); break;}      /* --injecttemplate              */
      case 'm': {strcpy(CLimportancedraws, optarg); break;}     /* --importanceresample          */
      case 'O': {strcpy(CLpriorpar, optarg); break;}            /* --priorparameters             */
      case 'h': {printhelpmessage(); return 0; break;}          /* --help                        */
    }
  }
  /* command line arguments require further processing...: */

  /* first of all figure out number of data sets / network size    */
  /* (taken to be size of either 'cachefile' or 'network' vector). */
  /* NOTE that much below hinges on this number.                   */
  if (CLcachefile[0] != '\0') {
    if (noVectorStrg(CLcachefile)) *coherentN = 1;
    else {
      parseCharacterOptionString(CLcachefile, &argNames, &argN);
      *coherentN = argN;
      for (i=0; i<argN; ++i) free(argNames[i]);
      free(argNames); 
    }
  }
  else if (CLnetwork[0] != '\0') {
    if (noVectorStrg(CLnetwork)) *coherentN = 1;
    else {
      parseCharacterOptionString(CLnetwork, &argNames, &argN);
      *coherentN = argN;
      for (i=0; i<argN; ++i) free(argNames[i]);
      free(argNames); 
    }
  }

  if (verbose) printf(" | Number of detectors: %d.\n", *coherentN);

  /* little sanity check: */
  if (! (((CLcachefile[0]!='\0') & (CLsimulateData[0]=='0')) 
         || ((CLnetwork[0]!='\0') & (CLsimulateData[0]=='1'))))
    printf(" : WARNING: strange parameter settings: either provide cache file(s) or have data simulated!\n");

  /* now allocate memory to hold data details: */
  DF = (DataFramework*) malloc(sizeof(DataFramework) * (*coherentN));

  /* initialise "DataFramework" details: */
  for (i=0; i<(*coherentN); ++i) {
    /*-- set dummy/default values for required command line arguments: -- */
    DF[i].data   = NULL;
    DF[i].window = NULL;
    DF[i].FTin   = NULL;
    DF[i].FTout  = NULL;
    DF[i].PsdEstimateStart = 0.0;       /* sec. */
    DF[i].PsdEstimateEnd   = 0.0;       /* sec. */
    DF[i].frameChannel = (char*) malloc(sizeof(char)*256);
    strcpy(DF[i].frameChannel, "\0");
    DF[i].datacachefile  = (char*) malloc(sizeof(char)*256);
    strcpy(DF[i].datacachefile, "\0");
    DF[i].noisecachefile = (char*) malloc(sizeof(char)*256);
    strcpy(DF[i].noisecachefile, "\0");
    DF[i].noiseModel = (char*) malloc(sizeof(char)*64);
    strcpy(DF[i].noiseModel, "initialLigo");
    DF[i].ifo = NULL;
  }

  tmpInt = (CLsimulateData[0]=='1');
  for (i=0; i<(*coherentN); ++i) DF[i].simulateData = tmpInt;
  /*--  fill in cache file names ("--cachefile" argument):  --*/
  if (CLcachefile[0] != '\0') {
    if (noVectorStrg(CLcachefile)) {
      if (*coherentN == 1) {
        strcpy(DF[i].datacachefile, CLcachefile);
        strcpy(DF[i].noisecachefile, CLcachefile);
      }
      else
        printf(" : WARNING: cannot make sense out of \"--cachefile\" argument!\n");
    }
    else {
      parseCharacterOptionString(CLcachefile, &argNames, &argN);      
      for (i=0; i<(*coherentN); ++i) {
        strcpy(DF[i].datacachefile, argNames[argN==1 ? 0 : i]);
        strcpy(DF[i].noisecachefile, argNames[argN==1 ? 0 : i]);
        /* (internally, data- and noise- cache files are still kept seperately.) */
      }
      for (i=0; i<argN; ++i) free(argNames[i]);
      free(argNames); 
    }
    if (verbose) {
      printf(" | cache files:\n");
      for (i=0; i<(*coherentN); ++i)
        printf(" | %d) \"%s\"\n", i+1, DF[i].datacachefile);
    }
  }
  else if (CLsimulateData[0]=='0') printf(" : WARNING: \"--cachefile\" argument missing!\n");

  /*--  "--network" argument:  --*/
  if (CLnetwork[0] != '\0') {
    if (noVectorStrg(CLnetwork)) {
      if (*coherentN == 1) {
        expandLocationString(CLnetwork, tmpChar1);
	DF[0].ifo = ifoPointer(ifolist, ifoN, tmpChar1);
      }
      else
        printf(" : WARNING: cannot make sense out of \"--network\" argument!\n");
    }
    else {
      parseCharacterOptionString(CLnetwork, &argNames, &argN);
      for (i=0; i<(*coherentN); ++i) {
        /* convert & copy location index to each DF */
        expandLocationString(argNames[argN==1 ? 0 : i], tmpChar1);
	DF[i].ifo = ifoPointer(ifolist, ifoN, tmpChar1);
      }
      for (i=0; i<argN; ++i) free(argNames[i]);
      free(argNames); 
    }
    if (verbose) {
      for (i=0; i<(*coherentN); ++i)
        printf(" | %d) \"%s\"\n", i+1, DF[i].ifo->name);
    }
  }

  /*--  "--tcenter" argument:  --*/
  if (CLtimeCenter[0] != '\0') {
    if (noVectorStrg(CLtimeCenter)) {
      tmpDbl = atof(CLtimeCenter);
      for (i=0; i<(*coherentN); ++i) DF[i].timeCenter = tmpDbl;
    }
    else {
      parseParameterOptionString(CLtimeCenter, &argNames, &argValues, &argN);
      if ((argN != *coherentN) && (argN != 1))
        printf(" : ERROR: incompatible number of \"--tcenter\" arguments (%d)!\n",argN);
      else for (i=0; i<(*coherentN); ++i) DF[i].timeCenter = argValues[argN==1 ? 0 : i];
      for (i=0; i<argN; ++i) free(argNames[i]);
      free(argNames); free(argValues);
    }
  }
  else printf(" : WARNING: \"--tcenter\" argument missing!\n");

  /*--  "--tbefore" argument:  --*/
  if (noVectorStrg(CLtimeBefore)) {
    tmpDbl = atof(CLtimeBefore);
    for (i=0; i<(*coherentN); ++i) DF[i].timeBefore = tmpDbl;
  }
  else {
    parseParameterOptionString(CLtimeBefore, &argNames, &argValues, &argN);
    if ((argN != *coherentN) && (argN != 1))
      printf(" : ERROR: incompatible number of \"--tbefore\" arguments (%d)!\n",argN);
    else for (i=0; i<(*coherentN); ++i) DF[i].timeBefore = argValues[argN==1 ? 0 : i];
    for (i=0; i<argN; ++i) free(argNames[i]);
    free(argNames); free(argValues);
  }

  /*--  "--tafter" argument:  --*/
  if (noVectorStrg(CLtimeAfter)) {
    tmpDbl = atof(CLtimeAfter);
    for (i=0; i<(*coherentN); ++i) DF[i].timeAfter = tmpDbl;
  }
  else {
    parseParameterOptionString(CLtimeAfter, &argNames, &argValues, &argN);
    if ((argN != *coherentN) && (argN != 1))
      printf(" : ERROR: incompatible number of \"--tafter\" arguments (%d)!\n",argN);
    else for (i=0; i<(*coherentN); ++i) DF[i].timeAfter = argValues[argN==1 ? 0 : i];
    for (i=0; i<argN; ++i) free(argNames[i]);
    free(argNames); free(argValues);
  }

  /*--  "--tukey" argument:  --*/
  tmpDbl = atof(CLtukeypar);
  for (i=0; i<(*coherentN); ++i) DF[i].tukeypar = tmpDbl;
 
  /*--  "--filechannel" argument:  --*/
  if (CLsimulateData[0]=='0') {
    if (CLframeChannel[0] != '\0') {
      if (noVectorStrg(CLframeChannel)) {
        for (i=0; i<(*coherentN); ++i) strcpy(DF[i].frameChannel, CLframeChannel);
        if (*coherentN > 1)
          printf(" : WARNING: same channel name (\"%s\") for all (%d) files. Correct?\n",
                 CLframeChannel, *coherentN);
      }
      else {
        parseCharacterOptionString(CLframeChannel, &argNames, &argN);
        if ((argN != *coherentN) && (argN != 1))
          printf(" : ERROR: incompatible number of \"--filechannel\" arguments (%d)!\n",argN);
        else for (i=0; i<(*coherentN); ++i) strcpy(DF[i].frameChannel, argNames[argN==1 ? 0 : i]);
        for (i=0; i<argN; ++i) free(argNames[i]);
        free(argNames);
      }
    }
    else printf(" : WARNING: \"--filechannel\" argument missing!\n");
  }

  /*--  "--psdestimatestart" argument:  --*/
  if (CLsimulateData[0]=='0') {
    if (CLPsdEstimateStart[0] != '\0') {
      if (noVectorStrg(CLPsdEstimateStart)) {
        tmpDbl = atof(CLPsdEstimateStart);
        for (i=0; i<(*coherentN); ++i) DF[i].PsdEstimateStart = tmpDbl;
      }
      else {
        parseParameterOptionString(CLPsdEstimateStart, &argNames, &argValues, &argN);
        if ((argN != *coherentN) && (argN != 1))
          printf(" : ERROR: incompatible number of \"--psdestimatestart\" arguments (%d)!\n",argN);
        else for (i=0; i<(*coherentN); ++i) DF[i].PsdEstimateStart = argValues[argN==1 ? 0 : i];
        for (i=0; i<argN; ++i) free(argNames[i]);
        free(argNames); free(argValues);
      }
    }
    else printf(" : WARNING: \"--psdestimatestart\" argument missing!\n");
  }

  /*--  "--psdestimateend" argument:  --*/
  if (CLsimulateData[0]=='0') {
    if (CLPsdEstimateStart[0] != '\0') {
      if (noVectorStrg(CLPsdEstimateEnd)) {
        tmpDbl = atof(CLPsdEstimateEnd);
        for (i=0; i<(*coherentN); ++i) DF[i].PsdEstimateEnd = tmpDbl;
      }
      else {
        parseParameterOptionString(CLPsdEstimateEnd, &argNames, &argValues, &argN);
        if ((argN != *coherentN) && (argN != 1))
          printf(" : ERROR: incompatible number of \"--psdestimateend\" arguments (%d)!\n",argN);
        else for (i=0; i<(*coherentN); ++i) DF[i].PsdEstimateEnd = argValues[argN==1 ? 0 : i];
        for (i=0; i<argN; ++i) free(argNames[i]);
        free(argNames); free(argValues);
      }
    }
    else printf(" : WARNING: \"--psdestimateend\" argument missing!\n");
  }

  /*--  "--fixspecdens" argument:  --*/
  tmpInt = (CLsimulatePsd[0]=='1');
  for (i=0; i<(*coherentN); ++i) DF[i].simulatePsd = tmpInt;

  /*--  "--psdestimaten" argument:  --*/
  if ((CLsimulateData[0]=='0') | (CLsimulatePsd[0]=='1')) {
    if (CLPsdEstimateN[0] != '\0') {
      if (noVectorStrg(CLPsdEstimateN)) {
        tmpInt = atoi(CLPsdEstimateN);
        for (i=0; i<(*coherentN); ++i) DF[i].PsdEstimateN = tmpInt;
      }
      else {
        parseParameterOptionString(CLPsdEstimateN, &argNames, &argValues, &argN);
        if ((argN != *coherentN) && (argN != 1))
          printf(" : ERROR: incompatible number of \"--psdestimaten\" arguments (%d)!\n",argN);
        else for (i=0; i<(*coherentN); ++i) DF[i].PsdEstimateN = ((int)argValues[argN==1 ? 0 : i]);
        for (i=0; i<argN; ++i) free(argNames[i]);
        free(argNames); free(argValues);
      }
    }
    else printf(" : WARNING: \"--psdestimaten\" argument missing!\n");
  }

  /*--  "--specdens" argument:  --*/
  if (CLsimulatePsd[0]=='1') {
    if (CLnoiseModel[0] != '\0') {
      if (noVectorStrg(CLnoiseModel)) {  /* (same noise model for all)        */
        for (i=0; i<(*coherentN); ++i) strcpy(DF[i].noiseModel, CLnoiseModel);
      }
      else {                             /* (individual noise model for each) */
        parseCharacterOptionString(CLnoiseModel, &argNames, &argN);
        if ((argN != *coherentN) && (argN != 1))
          printf(" : ERROR: incompatible number of \"--specdens\" arguments (%d)!\n",argN);
        else 
          for (i=0; i<(*coherentN); ++i) {
            strcpy(DF[i].noiseModel, argNames[argN==1 ? 0 : i]);
            if ((strcmp(argNames[i], "initialLigo") != 0) 
                & (strcmp(argNames[i], "advancedLigo") != 0)
                & (strcmp(argNames[i], "Virgo") != 0)
                & (strcmp(argNames[i], "Geo") != 0))
              printf(" : WARNING: \"--specdens\" option must be one of 'initialLigo', 'advancedLigo' or 'Virgo'!\n");
          }
        for (i=0; i<argN; ++i) free(argNames[i]);
        free(argNames); 
      }
    }
    else printf(" : WARNING: \"--specdens\" argument missing!\n");
  }

  /*--  "--freqlower" argument:  --*/
  if (noVectorStrg(CLfrequencyLower)) {
    tmpDbl = atof(CLfrequencyLower);
    for (i=0; i<(*coherentN); ++i) DF[i].minF = tmpDbl;
  }
  else {
    parseParameterOptionString(CLfrequencyLower, &argNames, &argValues, &argN);
    if ((argN != *coherentN) && (argN != 1))
      printf(" : ERROR: incompatible number of \"--freqlower\" arguments (%d)!\n",argN);
    else for (i=0; i<(*coherentN); ++i) DF[i].minF = argValues[argN==1 ? 0 : i];
    for (i=0; i<argN; ++i) free(argNames[i]);
    free(argNames); free(argValues);
  }

  /*--  "--frequpper" argument:  --*/
  if (noVectorStrg(CLfrequencyUpper)) {
    tmpDbl = atof(CLfrequencyUpper);
    for (i=0; i<(*coherentN); ++i) DF[i].maxF = tmpDbl;
  }
  else {
    parseParameterOptionString(CLfrequencyUpper, &argNames, &argValues, &argN);
    if ((argN != *coherentN) && (argN != 1))
      printf(" : ERROR: incompatible number of \"--frequpper\" arguments (%d)!\n",argN);
    else for (i=0; i<(*coherentN); ++i) DF[i].maxF = argValues[argN==1 ? 0 : i];
    for (i=0; i<argN; ++i) free(argNames[i]);
    free(argNames); free(argValues);
  }

  /*--  "--importanceresample" argument:  --*/
  importancedraws = atoi(CLimportancedraws);

  /*--  "--randomseed" argument:  --*/
  if (CLrandomseed[0] == '\0') printf(" : WARNING: no random seed specified!\n");
  else {
    if (noVectorStrg(CLrandomseed)) {
      seed1 = seed2 = atoi(CLrandomseed);
    }
    else {
      parseParameterOptionString(CLrandomseed, &argNames, &argValues, &argN);
      if (argN==1) {
        seed1 = seed2 = argValues[0];
      }
      else if (argN==2) {
        seed1 = argValues[0];
        seed2 = argValues[1];
      }
      else printf(" : WARNING: failed parsing 'randomseed' argument!\n");
      for (i=0; i<argN; ++i) free(argNames[i]);
      free(argNames); free(argValues);
    }
  }
  gsl_rng_set(GSLrandom, seed1);

  /*--  read / generate data:  --*/
  
  if (CLsimulateData[0]=='1') { /* generate fake data w/ given noise PSD:             */
    for (i=0; i<(*coherentN); ++i) {
      if (verbose) printf(" | Ifo %d: generating data using '%s' noise PSD.\n", i+1, DF[i].noiseModel);
      simulateData(&DF[i]);
    }
  }
  else {                        /* read data from files (and downsample by factor 4): */
    if (verbose) printf(" | reading data from files...\n");
    for (i=0; i<(*coherentN); ++i) {
      if (DF[i].datacachefile[0]!='\0'){ /* cache file supplied --> get file names from there */
        DF[i].datafilenames = cache2string(DF[i].datacachefile, &startGPS, &endGPS, tmpChar1);
        expandLocationString(tmpChar1, tmpChar2);
	DF[i].ifo = ifoPointer(ifolist, ifoN, tmpChar2);
        if ((startGPS > (DF[i].timeCenter - DF[i].timeBefore)) 
            | (endGPS < (DF[i].timeCenter + DF[i].timeAfter)))
          printf(" : ERROR: incomplete overlap between requested data and supplied GWF files!\n");
      }
      else printf(" : ERROR: need to supply (data) cache files!!\n");
    }
    for (i=0; i<(*coherentN); ++i)
      if (InitialisationOK) InitialisationOK = readData(&DF[i]);
  }

  /* set up window: */
  if (verbose) printf(" | using Tukey window with a=%.3f for all FTs.\n", DF[0].tukeypar);
  for (i=0; i<(*coherentN); ++i) {
    DF[i].window = (double*) malloc(sizeof(double)*DF[i].dataSize);
    DF[i].winss = 0.0;
    for (j=0; j<DF[i].dataSize; ++j) {
      DF[i].window[j] = tukeywindow(j, DF[i].dataSize, DF[i].tukeypar);
      DF[i].winss += DF[i].window[j] * DF[i].window[j];
    }
    if (verbose) printf(" | winss=%.1f, N=%.1f, winss/N=%.3f\n", 
                        DF[i].winss, ((double)DF[i].dataSize), DF[i].winss/((double)DF[i].dataSize));
  }

  /* initialize FT plan etc: */
  for (i=0; i<(*coherentN); ++i) {
    DF[i].FTin   = (double*) fftw_malloc(sizeof(double) * DF[i].dataSize);
    DF[i].FTSize = ((long) floor((((double)DF[i].dataSize)/2.0)+1.0));
    DF[i].FTout  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * DF[i].FTSize);
    DF[i].FTplan = fftw_plan_dft_r2c_1d(DF[i].dataSize, DF[i].FTin, DF[i].FTout, FFTW_ESTIMATE);
    DF[i].FTDeltaF = 1.0 / (DF[i].dataDeltaT * ((double)DF[i].dataSize));
    /* default settings for maxF (if still set to zero): */
    if (DF[i].maxF == 0.0) DF[i].maxF = (7.0/8.0) * ((DF[i].FTSize-1)*DF[i].FTDeltaF);
    /* = (7/8) x Nyquist frequency.  */
    /* sanity check (maxF < Nyquist ?): */
    if (DF[i].maxF > ((DF[i].FTSize-1)*DF[i].FTDeltaF))
      printf(" : ERROR: upper frequency bound above Nyquist frequency!!\n :        (%.1f > %.1f)\n",
	     DF[i].maxF, ((DF[i].FTSize-1)*DF[i].FTDeltaF));
    DF[i].minInd = ceil((DF[i].minF/((DF[i].FTSize-1)*DF[i].FTDeltaF))*(DF[i].FTSize-1));
    DF[i].maxInd = floor((DF[i].maxF/((DF[i].FTSize-1)*DF[i].FTDeltaF))*(DF[i].FTSize-1));
    if (verbose) {
      printf(" | Ifo %d (%s):\n", i+1, DF[i].ifo->name);
      printf(" | (time domain) data size: %d\n", DF[i].dataSize);
      printf(" | (freq. domain) FT size : %d\n", DF[i].FTSize);
      printf(" | Nyquist frequency      : %.1f Hz\n", (DF[i].FTSize-1)*DF[i].FTDeltaF);
      printf(" | frequency range: %.1f -- %.1f Hz\n", 
             DF[i].minInd*DF[i].FTDeltaF, DF[i].maxInd*DF[i].FTDeltaF);
    }
    /* Note that all Fourier transforms in the following also apply windowing to their input. */
    /* See also the "FTexec()" function.                                                      */
  
    /* Fourier transform (and also window) data: */
    DF[i].dataFT = (double complex*) malloc(sizeof(double complex)*DF[i].FTSize);
    FTexec(&DF[i], DF[i].data, DF[i].dataFT);
  }

  for (i=0; i<(*coherentN); ++i) {
    /* initialize noise power spectrum: */
    DF[i].powspec = (double*) malloc(sizeof(double) * DF[i].FTSize);

    /* fill in values for PSD... either estimate or use "known" spectra: */
    if (DF[i].simulateData) {
      if (DF[i].simulatePsd) {
        simulatePsdEstimation(&DF[i]);
      }
      else {
        if (verbose) printf(" | plugging in known spectrum for PSD\n");
        if (strcmp(DF[i].noiseModel, "advancedLigo")==0)
          for (j=0; j<DF[i].FTSize; ++j) 
            DF[i].powspec[j] = spectrumAdvLigo(((double)j)*DF[i].FTDeltaF);
        else if (strcmp(DF[i].noiseModel, "initialLigo")==0)
          for (j=0; j<DF[i].FTSize; ++j) 
            DF[i].powspec[j] = spectrumIniLigo(((double)j)*DF[i].FTDeltaF);
        else if (strcmp(DF[i].noiseModel, "Virgo")==0)
          for (j=0; j<DF[i].FTSize; ++j) 
            DF[i].powspec[j] = spectrumVirgo(((double)j)*DF[i].FTDeltaF);
        else 
          for (j=0; j<DF[i].FTSize; ++j) 
            DF[i].powspec[j] = spectrumGeo(((double)j)*DF[i].FTDeltaF);
      }
    }
    else {
      if (verbose) printf(" | estimating noise PSD from data\n");
      if (InitialisationOK) {
        if (DF[i].noisecachefile[0]!='\0'){ /* cache file supplied --> get file names from there */
          DF[i].noisefilenames = cache2string(DF[i].noisecachefile, &startGPS, &endGPS, tmpChar1);
          if (DF[i].PsdEstimateStart==0.0) DF[i].PsdEstimateStart = startGPS;
          if (DF[i].PsdEstimateEnd==0.0)   DF[i].PsdEstimateEnd   = endGPS;
          if ((startGPS>DF[i].PsdEstimateStart) | (endGPS<DF[i].PsdEstimateEnd))
            printf(" : ERROR: incomplete overlap between requested noise and supplied (noise) GWF files!\n");
        }
        else printf(" : ERROR: need to supply (noise) cache files!!\n");      
        InitialisationOK = estimatePSD(&DF[i]);
      }
    }
  }


  /*--  set up "MF" stuff:  --*/

  MF = (McmcFramework*) malloc(sizeof(McmcFramework));
  MF->studentDF = 3.0; /* (3 degrees of freedom: fixed setting for now) */

  /*--  "--iterations" argument:  --*/
  MF->iterations       = 1000000;
  if (CLiterations[0]!='\0') MF->iterations = atoi(CLiterations);

  /*--  "--template" argument:  --*/
  MF->template = char2template(CLtemplate);

  /*--  "--injecttemplate" argument:  --*/
  if (CLinjecttemplate[0] != '\0') 
    injecttemplate=char2template(CLinjecttemplate);
  else 
    injecttemplate = MF->template;
  /* (by default injection template same as recovery template) */

  MF->logfilename = (char*) malloc(sizeof(char)*256);
  if (CLlogfilename[0] != '\0')
    strcpy(MF->logfilename, CLlogfilename);
  else printf(" : WARNING: \"--logfile\" argument missing!\n");

  /* derive GMST:                                                                  */
  /* first average over individual `timeCenter' values (usually identical anyway): */
  tmpDbl = 0.0;
  for (i=0; i<*coherentN; ++i)
    tmpDbl += DF[i].timeCenter - DF[0].timeCenter;
  tmpDbl /= *coherentN;
  tmpDbl += DF[0].timeCenter;
  MF->GMST = GMST(tmpDbl);
  if (verbose) printf(" | reference GMST is %.4f (rad).\n",MF->GMST);

  MF->celestialCoords = 1; /* log right ascension & declination (instead of latitude & longitude) */

  /* set corresponding parameter set etc: */

  /*--- Inspiral, no spin, 9 parameters: ---*/
  if ((MF->template == iR2PN)
      | (MF->template == i20SP)
      | (MF->template == i25SP)
      | (MF->template == i2025)
      | (MF->template == i2535)
      | (MF->template == iLALTT2PN00)
      | (MF->template == iLALTT2PN10)
      | (MF->template == iLALTT2PN15)
      | (MF->template == iLALTT2PN20)
      | (MF->template == iLALTT3PN00)
      | (MF->template == iLALTT3PN10)
      | (MF->template == iLALTT3PN15)
      | (MF->template == iLALTT3PN20)
      | (MF->template == iLALIMRPhenomA)
      | (MF->template == iLALEOBNR)     ) {
    MF->parameterset=InspiralNoSpin;
    /* (already) set up starting value vector:                   */
    /* (proper list of variable names (dimension!) required for  */
    /* proposal covariance matrix, fixed values, &c. below)      */
    vectorInit(&MF->startvalue);
    vectorSetup(&MF->startvalue, MF->parameterset);
    MF->pardim = MF->startvalue.dimension;
    if (verbose) 
      printf(" | setting up MCMC for template '#%d' / parameterset '#%d', %d parameters.\n",
             MF->template, MF->parameterset, MF->pardim);
    /* set up covariance matrix: */
    MF->covariance = (double*) malloc(sizeof(double) * (MF->pardim*MF->pardim));
    for (i=0; i<(MF->pardim*MF->pardim); ++i) MF->covariance[i] = 0.0;
    /* set variances (via standard deviations): */
    setstdev(MF, parindex(MF,"latitude"),     1.5 * pi/180.0); /* <-- relevant for sky location! */
    setstdev(MF, parindex(MF,"longitude"),    1.0 * pi/180.0); /* <-- this value is ignored!     */
    setstdev(MF, parindex(MF,"polarisation"), 5.0 * pi/180.0); /*     (see proposal function)    */
    setstdev(MF, parindex(MF,"inclination"),  5.0 * pi/180.0);
    setstdev(MF, parindex(MF,"phase"),        20.0 * pi/180.0);
    setstdev(MF, parindex(MF,"time"),         0.001);
    setstdev(MF, parindex(MF,"logdistance"),  0.1);
    setstdev(MF, parindex(MF,"chirpmass"),    /*0.002);*/ 0.0001);
    setstdev(MF, parindex(MF,"massratio"),    0.005);
    /* set covariances (via correlations): */
    setcor(MF, parindex(MF,"chirpmass"),   parindex(MF,"massratio"),   0.90);
    setcor(MF, parindex(MF,"massratio"),   parindex(MF,"time"),        0.80);
    setcor(MF, parindex(MF,"logdistance"), parindex(MF,"inclination"), 0.70);
    setcor(MF, parindex(MF,"chirpmass"),   parindex(MF,"time"),        0.60);
    /* set (6) prior parameters: */
    vectorInit(&MF->priorparameters);
    /* default settings:         */
    vectorAdd(&MF->priorparameters, "massLower", 1.0);
    vectorAdd(&MF->priorparameters, "massUpper", 15.0);
    vectorAdd(&MF->priorparameters, "timeLower", DF[0].timeCenter-0.050);
    vectorAdd(&MF->priorparameters, "timeUpper", DF[0].timeCenter+0.050);
    vectorAdd(&MF->priorparameters, "dist90", 40.0);
    vectorAdd(&MF->priorparameters, "dist10", 80.0);
    if (CLpriorpar[0] != '\0') {
      parseParameterOptionString(CLpriorpar, &argNames, &argValues, &argN);
      if (argN == 6)
        for (i=0; i<6; ++i)
          MF->priorparameters.value[i] = argValues[i];
      else
        printf(" : ERROR: incorrect number (%d) of prior parameters specified!", argN);
      for (i=0; i<argN; ++i) free(argNames[i]);
      free(argNames); free(argValues);
    }
    if (verbose){
      printf(" | prior settings:\n");
      printf(" |   mass range          : %.1f -- %.1f\n", 
             vectorGetValue(&MF->priorparameters,"massLower"), 
             vectorGetValue(&MF->priorparameters,"massUpper"));
      printf(" |   time range          : %.3f -- %.3f\n", 
             vectorGetValue(&MF->priorparameters,"timeLower"),
             vectorGetValue(&MF->priorparameters,"timeUpper"));
      printf(" |   2-2-Ms 90%% detection: %.1f Mpc\n", 
             vectorGetValue(&MF->priorparameters,"dist90"));
      printf(" |   2-2-Ms 10%% detection: %.1f Mpc\n", 
             vectorGetValue(&MF->priorparameters,"dist10"));
    }    
  }

  /*--- sine-gaussian burst, 8 parameters: ---*/
  else if (MF->template == bSineGaussian) {
    MF->parameterset = BurstSineGaussian;
    /* (already) set up starting value vector:                   */
    /* (proper list of variable names (dimension!) required for  */
    /* proposal covariance matrix, fixed values, &c. below)      */
    vectorInit(&MF->startvalue);
    vectorSetup(&MF->startvalue, MF->parameterset);
    MF->pardim = MF->startvalue.dimension;
    /* set up covariance matrix: */
    MF->covariance = (double*) malloc(sizeof(double) * (MF->pardim*MF->pardim));
    for (i=0; i<(MF->pardim*MF->pardim); ++i) MF->covariance[i] = 0.0;
    /* set variances: */
    setstdev(MF, parindex(MF,"latitude"),      1.0 * pi/180.0);
    setstdev(MF, parindex(MF,"longitude"),     1.0 * pi/180.0);
    setstdev(MF, parindex(MF,"polarisation"),  1.0 * pi/180.0);
    setstdev(MF, parindex(MF,"phase"),         5.0 * pi/180.0);
    setstdev(MF, parindex(MF,"time"),          0.001);
    setstdev(MF, parindex(MF,"logamplitude"),  0.2);
    setstdev(MF, parindex(MF,"logsigma"),      0.1);
    setstdev(MF, parindex(MF,"frequency"),     0.1);
    /* set covariances: */
    setcor(MF, parindex(MF,"phase"), parindex(MF,"time"), -0.9);
    /* set (6) prior parameters: */
    vectorInit(&MF->priorparameters);
    /* default settings: */
    vectorAdd(&MF->priorparameters, "freqLower",   1.0);
    vectorAdd(&MF->priorparameters, "freqUpper",   1500.0);
    vectorAdd(&MF->priorparameters, "timeLower",   DF[0].timeCenter-0.050);
    vectorAdd(&MF->priorparameters, "timeUpper",   DF[0].timeCenter+0.050);
    vectorAdd(&MF->priorparameters, "ampliExpect", 1.0e-20);
    vectorAdd(&MF->priorparameters, "sigmaExpect", 0.050);
    if (CLpriorpar[0] != '\0') {
      parseParameterOptionString(CLpriorpar, &argNames, &argValues, &argN);
      if (argN == 6)
        for (i=0; i<6; ++i)
          MF->priorparameters.value[i] = argValues[i];
      else
        printf(" : ERROR: incorrect number (%d) of prior parameters specified!", argN);
      for (i=0; i<argN; ++i) free(argNames[i]);
      free(argNames); free(argValues);
    }

    if (verbose){
      printf(" | prior settings:\n");
      printf(" |   frequency range    : %.1f -- %.1f\n", 
             vectorGetValue(&MF->priorparameters,"freqLower"), 
             vectorGetValue(&MF->priorparameters,"freqUpper"));
      printf(" |   time range         : %.3f -- %.3f\n", 
             vectorGetValue(&MF->priorparameters,"timeLower"), 
             vectorGetValue(&MF->priorparameters,"timeUpper"));
      printf(" |   expected amplitude : %.3e\n", 
             vectorGetValue(&MF->priorparameters,"ampliExpect"));
      printf(" |   expected sigma     : %.3f\n", 
             vectorGetValue(&MF->priorparameters,"sigmaExpect"));
    }
  }

  else
    printf(" : ERROR: referring to undefined template waveform (#%d) in 'initMF()'!\n",
           MF->template);

  /* determine Cholesky decomposition of cocariance matrix for proposal functions: */
  MF->covcholesky = gsl_matrix_alloc(MF->pardim, MF->pardim);
  for (i=0; i<MF->pardim; ++i)
    for (j=0; j<MF->pardim; ++j)
      gsl_matrix_set(MF->covcholesky, i, j, MF->covariance[i+j*MF->pardim]);
  gsl_linalg_cholesky_decomp(MF->covcholesky);

  /* determine to be fixed & to be inferred parameters:  */
  vectorInit(&MF->fixed);
  if (CLfixedpar[0] != '\0') {
    parseParameterOptionString(CLfixedpar, &argNames, &argValues, &argN);
    for (i=0; i<argN; ++i) {
      if (vectorIsElement(&MF->startvalue, argNames[i]))
        vectorAdd(&MF->fixed, argNames[i], argValues[i]);
      else if (strcmp(argNames[i],"declination")==0)
        vectorAdd(&MF->fixed, "latitude", argValues[i]);
      else if (strcmp(argNames[i],"rightascension")==0)
        vectorAdd(&MF->fixed, "longitude", longitude(argValues[i],MF->GMST));
      else {
        printf(" : WARNING: unknown parameter \"%s\" provided as value to be fixed!\n", argNames[i]); 
        printf(" :          (allowed parameters:");
        printf(" \"declination\", \"rightascension\",");
        for (j=0; j<MF->startvalue.dimension; ++j) {
          printf(" \"%s\"", MF->startvalue.name[j]);
          if (j<(MF->startvalue.dimension-1)) printf(",");
          else printf(")\n");
	}
      }
    }
    for (i=0; i<argN; ++i) free(argNames[i]);
    free(argNames); free(argValues);
    /* printf(" : fixed parameters:\n"); vectorPrint(&MF->fixed); */
  }

  /* determine starting values (if provided in "--start" argument): */
  /* ('startvalue' vector was initialised above)                    */
  priordraw(MF, &MF->startvalue);
  if (CLstartpar[0] != '\0') {
    parseParameterOptionString(CLstartpar, &argNames, &argValues, &argN);
    for (i=0; i<argN; ++i) {
      if (vectorIsElement(&MF->startvalue, argNames[i]))
        vectorSetValue(&MF->startvalue, argNames[i], argValues[i]);
      else if (strcmp(argNames[i],"declination")==0)
        vectorSetValue(&MF->startvalue, "latitude", argValues[i]);
      else if (strcmp(argNames[i],"rightascension")==0)
        vectorSetValue(&MF->startvalue, "longitude", longitude(argValues[i],MF->GMST));
      else {
        printf(" : WARNING: unknown parameter \"%s\" provided as starting value!\n", argNames[i]); 
        printf(" :          (allowed parameters:");
        printf(" \"declination\", \"rightascension\",");
        for (j=0; j<MF->startvalue.dimension; ++j) {
          printf(" \"%s\"", MF->startvalue.name[j]);
          if (j<(MF->startvalue.dimension-1)) printf(",");
          else printf(")\n");
	}
      }
    }
    for (i=0; i<argN; ++i) free(argNames[i]);
    free(argNames); free(argValues);
    if (importancedraws > 0)
      printf(" : WARNING: provided starting values are going to be ignored!\n");
  }

  if (CLguesspar[0] != '\0') {
    parseParameterOptionString(CLguesspar, &argNames, &argValues, &argN);
    MF->guessparameters = (double*) malloc(sizeof(double)*argN);
    for (j=0; j<argN; ++j) {
      MF->guessparameters[j] = argValues[j];
    }
    for (i=0; i<argN; ++i) free(argNames[i]);
    free(argNames); free(argValues);
    if (verbose) {
      printf(" | %d parameter guesses / trigger values provided:\n", argN);
      printf(" |      chirpmass :  %.4f Msun\n", MF->guessparameters[0]);
      printf(" |      massratio :  %.4f\n", MF->guessparameters[1]);
      printf(" |      time      :  %.4f s\n", MF->guessparameters[2]);
      printf(" |      distance  :  %.2f Mpc\n", MF->guessparameters[3]);
    }
  }
  else MF->guessparameters = NULL;

  /* check for discrepancies between starting values and fixed values: */
  for (i=0; i<MF->fixed.dimension; ++i)
    if (vectorGetValue(&MF->startvalue,MF->fixed.name[i]) != MF->fixed.value[i]) {
      printf(" : WARNING: discrepancy between starting value and fixed value!\n");
      printf(" :          (\"%s\" is %.3e but should be %.3e)\n",
             MF->fixed.name[i], 
             vectorGetValue(&MF->startvalue,MF->fixed.name[i]),
             MF->fixed.value[i]);
    }

  /* determine injection parameters & inject signal: */
  vectorInit(&injectpar);
  if (CLinjectpar[0] != '\0') {
    vectorSetup(&injectpar, MF->parameterset);
    parseParameterOptionString(CLinjectpar, &argNames, &argValues, &argN);
    /* set parameter values; celestial coordinates are instantly translated to geographical. */
    for (i=0; i<argN; ++i) {
      if (vectorIsElement(&injectpar, argNames[i]))
        vectorSetValue(&injectpar, argNames[i], argValues[i]);
      else if (strcmp(argNames[i],"declination")==0)
        vectorSetValue(&injectpar, "latitude", argValues[i]);
      else if (strcmp(argNames[i],"rightascension")==0)
        vectorSetValue(&injectpar, "longitude", longitude(argValues[i],MF->GMST));
      else {
        printf(" : WARNING: unknown parameter \"%s\" provided as injection value!\n", argNames[i]); 
        printf(" :          (allowed parameters:");
        printf(" \"declination\", \"rightascension\",");
        for (j=0; j<injectpar.dimension; ++j) {
          printf(" \"%s\"", injectpar.name[j]);
          if (j<(injectpar.dimension-1)) printf(",");
          else printf(")\n");
	}
      }
    }
    for (i=0; i<argN; ++i) free(argNames[i]);
    free(argNames); free(argValues);
    if (verbose) {
      printf(" : injection:\n"); vectorPrint(&injectpar);
    }
    inject(DF, *coherentN, injecttemplate, &injectpar);
    /* dumptemplates(DF, &injectpar, 
                     "/home/christian/temp/templatesF.txt", 
                     "/home/christian/temp/templatesT.txt"); */
    vectorDispose(&injectpar);
  }

  MF->secPerIteration = 0.0;

  /* for (i=0; i<(*coherentN); ++i) printDF(&DF[i]); */

  long lhdf = 0;
  for (i=0; i<*coherentN; ++i)
    lhdf += 2 * (DF[i].maxInd - DF[i].minInd + 1);
  printf(" : (log-) likelihood degrees-of-freedom: %d\n",lhdf);

  /* initialise random number gernerator for (now following) MCMC part: */
  gsl_rng_set(GSLrandom, seed2);

  if (importancedraws > 0)
    importanceresample(DF, *coherentN, MF, &MF->startvalue, importancedraws, 1);

  *DFarg = DF;
  *MFarg = MF;
    
  return InitialisationOK;
} /*--  end of "init()".  --*/


void clearDF(DataFramework *DF)
/* cleans up pointers etc. in a "DataFramework" structure. */
{
  if (verbose) printf(" | cleaning up 'DataFramework' structure ...");
  free(DF->data);
  free(DF->window);
  fftw_free(DF->FTin);
  fftw_free(DF->FTout);
  fftw_destroy_plan(DF->FTplan);
  free(DF->powspec);
  free(DF->dataFT);
  free(DF);
  DF = NULL;
  if (verbose) printf(" done.\n");
}


int readData(DataFramework *DF)
/* this function reads the data, downsamples and windows it,    */
/* copies the result to "DF->data", and sets all implied values */
/* (DF->data, DF->dataStart, DF->dataDeltaT, DF->dataSize).     */
/* Requires frame file settings etc already to be given         */
/* in corresponding elements of "DF".                           */
{
  double from, to, extramargin;
  struct FrFile *iFile=NULL;
  struct FrVect *dataVector=NULL;
  double *origData=NULL; 
  int samplerate;
  int i, N;

  /* first determine the time range aimed for AFTER downsampling: */
  if (forceFlatTukey)
    extramargin = (DF->timeBefore+DF->timeAfter) * 0.5 * (DF->tukeypar/(1.0-DF->tukeypar));
  else
    extramargin = 0.0;
  from  = floor(DF->timeCenter - DF->timeBefore - extramargin);
  to    =  ceil(DF->timeCenter + DF->timeAfter  + extramargin);
  /* this range is rounded to an integer number of seconds */
  /* in order to make life easier for the DFT algorithm.   */

  if (downsamplingfactor > 1){
    /* downsampling will shorten the above range, so another extra margin is added. */
    /* This again requires the sampling rate to be known,                           */
    /* so a few samples are read in order to check:                                 */
    iFile = FrFileINew(DF->datafilenames);
    if (iFile == NULL) {
      printf(" : ERROR opening frame file!\n");
      return 0;
    }

    dataVector = FrFileIGetVectF(iFile, DF->frameChannel, floor(DF->timeCenter), 1.0);
    if (dataVector == NULL) {
      printf(" : ERROR reading frame file!\n");
      return 0;
    }
    samplerate = (1.0 / (dataVector->dx[0]) +0.5); /* samples per second                        */
    FrVectFree(dataVector);                        /* (add 0.5 for correct truncation/rounding) */
    FrFileIEnd(iFile);
    if (verbose) printf(" | data sampling rate: %d Hz\n", samplerate);
    /* finished determining sampling rate. */

    /* little sanity check: */
    if (DF->maxF > (0.11117*((double)samplerate)))
      printf(" :   !! WARNING !!\n :   'maxF' lies outside downsampling-filter passband !\n :   (%.1f > %.1f)\n", DF->maxF, (0.11117*((double)samplerate)));
    /* the figure "0.11117" refers to the hard-coded filter band settings */
    /* within the "filter()" function.                                    */

    /* account for how much the downsampling/filtering will nibble off the data */
    /* (downsamplingfactor*20, according to LAL documentation):                 */
    extramargin = (downsamplingfactor*20) / ((double) samplerate);
    from -= extramargin;
    to   += extramargin;
  }

  DF->rawDataRange = to-from; /* (save this figure for re-use for noise spectrum estimation) */
  if (verbose)
    printf(" | time range to be read: %.5f -- %.5f (%.5f sec.)\n", from, to, DF->rawDataRange);

  /*-- open frame file(s): --*/
  iFile = FrFileINew(DF->datafilenames);
  if (iFile == NULL) {
    printf(" : ERROR opening frame file(s)!\n");
    return 0;
  }

  /*-- read data: --*/
  /* if (doubleprecision)
     dataVector = FrFileIGetVectD(iFile, DF->frameChannel, from, DF->rawDataRange);
     else
     dataVector = FrFileIGetVectF(iFile, DF->frameChannel, from, DF->rawDataRange); */
  dataVector = FrFileIGetVectF(iFile, DF->frameChannel, from, DF->rawDataRange);
  if (dataVector == NULL) {
    printf(" : ERROR reading frame file(s)!\n");
    return 0;
  }

  N = dataVector->nData;  /* number of samples read from file(s)   */
  samplerate = (1.0 / (dataVector->dx[0]) +0.5); /* samples per second;                       */
                                                 /* (add 0.5 for correct truncation/rounding) */

  /*-- allocate memory for raw data: --*/
  origData = malloc(sizeof(double) * N);
  /*-- then copy over:               --*/
  for (i=0; i<N; ++i)
    origData[i] = dataVector->dataF[i];

  /*-- release the FrVect objects and close frame file(s) --*/
  FrVectFree(dataVector);
  FrFileIEnd(iFile);

  if (downsamplingfactor > 1){
    /*-- downsample (by factor 4)   --*/
    /*-- !! changes value of `N' !! --*/
    if (verbose) printf(" | downsampling... ");
    DF->data = downsample(origData, &N, downsamplingfactor);
    samplerate /= downsamplingfactor;
    if (verbose) {
      printf("new sample size and rate: %d at %d Hz\n", N, samplerate);
      /*printf(" | (double-check:  %d = %d x %d + %d)\n", N, N/samplerate, samplerate, N % samplerate); */
    }
    /* set parameters: */
    /*DF->dataStart  = from + ((double)(ncoef-1))/((double)(samplerate*4)); */
    DF->dataStart  = from + ((double)(downsamplingfactor*20))/((double)(samplerate*downsamplingfactor)); 
    /*  --> need to account for sampling rate BEFORE downsampling, hence factor 4  */
    DF->dataDeltaT = 1.0/((double) samplerate);
    DF->dataSize   = N;
  }
  else { /*-- no downsampling --*/
    if (verbose) printf(" | (no downsampling)\n");
    DF->data = (double*) malloc(sizeof(double) * N);
    for (i=0; i<N; ++i) DF->data[i] = origData[i];
    /* set parameters: */
    DF->dataStart  = from;
    DF->dataDeltaT = 1.0/((double) samplerate);
    DF->dataSize   = N;
  }

  free(origData);
  return 1;
}


char *cache2string(char *cachefilename, double *startTime, double *endTime, char *locationstrg)
/*****************************************************************/
/*  Reads a cache files and returns a Frame-library-digestible   */
/*  character string of GWF filenames.                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*  Expects cache file to look something like:                   */
/*    [...]                                                      */
/*    L L1_RDS_C03_L2 847546596 128                              */
/*    file://localhost/nfsdata/nfsdata24/S5/L1_RDS_C03_L2/L/847543000-847552999/L-L1_RDS_C03_L2-847546596-128.gwf */
/*    L L1_RDS_C03_L2 847546724 128                              */
/*    file://localhost/nfsdata/nfsdata24/S5/L1_RDS_C03_L2/L/847543000-847552999/L-L1_RDS_C03_L2-847546724-128.gwf */
/*    [...]                                                      */
/*  I.e.: two character strings, then start time and duration    */
/*  of following frame file, and eventually the actual filename  */
/*  (preceded by "file://localhost") for each entry.             */
/*  The returned character string contains the individual        */
/*  filenames, separated by white space,                         */
/*  as required by the Frame library's "FrFileINew()" function.  */
/*  "startTime" and "endTime" return the GPS time range          */
/*  spanned by the files.                                        */
/*  "locationstrg" is the first character(-string) in the first  */
/*  line of the file, and is used to identify the ifo.           */
/*****************************************************************/
{
  FILE *cachefile;
  int i, j=1;
  char location[8];
  char dataname[32];
  char errorstrg[17];
  int filenameerror = 0;
  int startGPS, deltaGPS;
  char filename[256];      /* (this makes 256 characters the maximal length of an individual file name. Feel free to change if necessary but keep track of implied changes below.) */
  char *allfilenames=NULL; /* vector of filenames  */
  char *cptr=NULL;         /* dummy pointer        */
  long nchar=0;            /* number of characters */
  /* initialise: */
  allfilenames = (char*) malloc(sizeof(char));
  allfilenames[0] = '\0';
  cachefile = fopen(cachefilename, "r");
  /* loop over cache file's entries: */
  j = fscanf(cachefile, "%s %s %d %d %s", 
             locationstrg, dataname, &startGPS, &deltaGPS, filename);
  *startTime = startGPS;
  while (j != EOF) {
    *endTime = startGPS + deltaGPS;
    /* check length of new filename: */
    i=0;
    while ((i<256) && (filename[i] != '\0')) ++i;
    if (filename[i] != '\0') printf(" : WARNING: encountered un-terminated file name string in 'cache2string()'.\n");
    /* (i gives number of 'non-zero' characters, or index of end-of-string character)       */
    /* issue error message if string is too short to even hold the "file://localhost" part, */
    /* or if that part is different:                                                        */
    if ((i <= 16) || (strncmp(filename,"file://localhost",16)!=0)) {
      strncpy(errorstrg, filename, 16);
      errorstrg[16] = '\0';
      printf(" : ERROR: invalid file name string in 'cache2string()'.\n");
      printf(" :        cache file :  \"%s\"\n", cachefilename);
      printf(" :        expected   :  \"file://localhost...\"\n");
      printf(" :        found      :  \"%s...\"\n", errorstrg);
      filenameerror = 1;
    }
    /* allocate new character vector:                */
    cptr = allfilenames;
    allfilenames = (char*) malloc(sizeof(char) * (nchar + (i-16) + (nchar!=0) + 1));
    /* copy old names over:                          */
    for (i=0; i<nchar; ++i) {
      allfilenames[i] = cptr[i];
    }
    /* insert white space between old and new names: */
    if (nchar!=0) allfilenames[nchar] = ' ';
    /* copy over new name (skip the "file:/" part):  */
    i = 16;
    while ((i<256) && (filename[i] != '\0')) {
      allfilenames[nchar + (nchar!=0) + (i-16)] = filename[i];
      ++i;
    }
    allfilenames[nchar+(nchar!=0)+(i-16)] = '\0';
    nchar = nchar + (nchar!=0) + (i-16);
    free(cptr);
    /* (attempt to) read next cache file entry:      */
    j = fscanf(cachefile, "%s %s %d %d %s", 
               location, dataname, &startGPS, &deltaGPS, filename);
  }
  if (filenameerror)
    printf(" :        resulting file name string: \"%s\"\n",allfilenames);
  fclose(cachefile);
  return allfilenames;
}


void generateNoise(double deltat, int N, 
                   double* noise, char *model)
/* "deltat" and "N" are resolution and length of requested noise vector. */
/* The spectral density (e.g. "spectrumIniLigo(f)" below) is the         */
/* one-sided (!) power spectral density of the resulting data.           */
{
  int Neven = (N % 2 == 0);
  int halfN = Neven ? N/2 : (N-1)/2;
  double deltaf = 1.0 / (((double)N) * deltat);
  int i;
  double f, kappa, stdev, a, b, real, imag;
  fftw_complex *FTinput=NULL;
  fftw_plan InvFTplan;
  int modelIndex = 0;

  FTinput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(halfN+1));

  if (strcmp(model, "initialLigo")==0)  
    modelIndex = 0; /* initial Ligo  */
  else if (strcmp(model, "advancedLigo")==0)  
    modelIndex = 1; /* advanced Ligo */
  else if (strcmp(model, "Virgo")==0)  
    modelIndex = 2; /* Virgo         */
  else if (strcmp(model, "Geo")==0)  
    modelIndex = 3; /* Geo           */
  else printf(" : ERROR: referring to undefined noise model ('%s') in 'generateNoise()'!\n",model);

  for (i=0; i<=halfN; ++i){
    f = ((double) i) * deltaf;
    kappa = (i==0) | (Neven & (i==halfN)) ? 0.0 : 1.0;
    if (modelIndex==0)
      stdev = exp(0.5*(spectrumIniLigo(f) - 2.0*log(1.0+kappa)));
    else if (modelIndex==1)
      stdev = exp(0.5*(spectrumAdvLigo(f) - 2.0*log(1.0+kappa)));
    else if (modelIndex==2)
      stdev = exp(0.5*(spectrumVirgo(f) - 2.0*log(1.0+kappa)));
    else if (modelIndex==3)
      stdev = exp(0.5*(spectrumGeo(f) - 2.0*log(1.0+kappa)));
    /*
    a = gennor(0.0, stdev);
    b = gennor(0.0, stdev) * kappa;
    */
    a = gsl_ran_ugaussian(GSLrandom) * stdev;
    b = gsl_ran_ugaussian(GSLrandom) * stdev * kappa;
    real = sqrt(((double)N)/deltat) * a;
    imag = -sqrt(((double)N)/deltat) * b;
    FTinput[i] = real + I*imag;
  }
  InvFTplan = fftw_plan_dft_c2r_1d(N, FTinput, noise, FFTW_ESTIMATE);
  fftw_execute(InvFTplan);
  fftw_destroy_plan(InvFTplan);
  fftw_free(FTinput);
  for (i=0; i<N; ++i)
    noise[i] /= N;
}


void simulateData(DataFramework *DF)
/* requires some settings from "DF", in particular: */
/* DF->timeCenter, DF->timeBefore, DF->timeAfter,   */
/* DF->dataDeltaT, DF->noiseModel                   */
{
  double from, to, extramargin;

  /* set sampling rate: */
  if (strcmp(DF->noiseModel, "initialLigo")==0)  
    DF->dataDeltaT = 1.0/16384.0; /* initial Ligo  */
  else if (strcmp(DF->noiseModel, "advancedLigo")==0)  
    DF->dataDeltaT = 1.0/16384.0; /* advanced Ligo */
  else  if (strcmp(DF->noiseModel, "Virgo")==0)  
    DF->dataDeltaT = 1.0/20000.0; /* Virgo         */
  else  if (strcmp(DF->noiseModel, "Geo")==0)  
    DF->dataDeltaT = 1.0/16384.0; /* Geo           */
  else printf(" : ERROR: referring to undefined noise model ('%s') in 'simulateData()'!\n",DF->noiseModel);
  DF->dataDeltaT *= downsamplingfactor;

  if (forceFlatTukey)
    extramargin = (DF->timeBefore+DF->timeAfter) * 0.5 * (DF->tukeypar/(1.0-DF->tukeypar));
  else
    extramargin = 0.0;
  from  = floor(DF->timeCenter - DF->timeBefore - extramargin);
  to    =  ceil(DF->timeCenter + DF->timeAfter  + extramargin);
  DF->dataStart = from;
  DF->dataSize = ((long) ((to-from) / DF->dataDeltaT));
  DF->data = (double*) malloc(sizeof(double)*DF->dataSize);
  generateNoise(DF->dataDeltaT, DF->dataSize, DF->data, DF->noiseModel);
}


int estimatePSD(DataFramework *DF)
/* estimates the 1-sided power spectral density (PSD) from the data; */
/* this function requires "readData()" to have run first,            */
/* as it relies on some of its settings...                           */
/* The FFT-related elements (FTin, FTout, FTplan,...) also need      */
/* to be initialised already.                                        */
{
  struct FrFile *iFile=NULL;
  struct FrVect *dataVector=NULL;
  int i, N, nsegment;
  double from, step;
  double *origData=NULL;
  double *dsdata=NULL;
  double complex *dataFT=NULL;
  double power, psdcoef;
  int PsdOk = 1;

  iFile = FrFileINew(DF->noisefilenames);
  if (iFile == NULL) {
    printf(" : ERROR opening noise file(s)!\n");
    PsdOk = 0;
  }

  from = DF->PsdEstimateStart;     /* starting point of 1st segment     */
  step = ceil(DF->rawDataRange);   /* stepwidth for subsequent segments */
  for (i=0; i<DF->FTSize; ++i) DF->powspec[i] = 0.0;

  /* loop over data segments: */
  nsegment = 0;
  while (((from + DF->rawDataRange) < DF->PsdEstimateEnd)
         & (nsegment < DF->PsdEstimateN)) {
    /* (this loop stops when either "DF->PsdEstimateEnd" is reached,          */
    /* or "DF->PsdEstimateN" segments have been read, whichever occurs first) */

    /* check for overlap with actual data: */
    if (((DF->timeCenter-DF->timeBefore > from) 
         & (DF->timeCenter-DF->timeBefore < from+DF->rawDataRange))
        | ((DF->timeCenter+DF->timeAfter > from) 
           & (DF->timeCenter+DF->timeAfter < from+DF->rawDataRange))
        | ((DF->timeCenter-DF->timeBefore < from) 
           & (DF->timeCenter+DF->timeAfter > from+DF->rawDataRange)))
      printf(" : WARNING: spectrum estimation data overlaps with data to be analysed !!\n");

    /*-- read data: --*/
    dataVector = FrFileIGetVectF(iFile, DF->frameChannel, from, DF->rawDataRange);
    if (dataVector == NULL) {
      printf(" : ERROR reading noise file(s)!\n");
      printf(" : (trying to read channel \"%s\", GPS %.3f--%.3f)\n", DF->frameChannel, from, from+DF->rawDataRange);
      printf(" : File list: %s\n", DF->noisefilenames);
      PsdOk = 0;
    }
    N = dataVector->nData;  /* number of samples read from file(s) */
    if (origData==NULL)
      origData = (double*) malloc(sizeof(double) * N);
    for (i=0; i<N; ++i)
      origData[i] = dataVector->dataF[i];
    FrVectFree(dataVector);
    if (downsamplingfactor > 1){
      /*-- downsample data: --*/
      dsdata = downsample(origData, &N, downsamplingfactor);
    }
    else {
      dsdata = (double*) malloc(sizeof(double) * N);
      for (i=0; i<N; ++i) dsdata[i] = origData[i];
    }
    /*-- window and FT data: --*/
    if (dataFT==NULL)
      dataFT = (double complex*) malloc(sizeof(double complex) * (N/2 + 1));
    FTexec(DF, dsdata, dataFT);
    for (i=0; i<DF->FTSize; ++i) {
      power = cabs(dataFT[i]);
      power *= power;
      DF->powspec[i] += power;
    }
    free(dsdata);
    from += step;
    nsegment += 1;
  }
  FrFileIEnd(iFile);
  free(origData);
  free(dataFT);

  psdcoef = log(2.0) + log(DF->dataDeltaT) - log(DF->dataSize) - log(nsegment);
  /* the factor  dataDeltaT/dataSize  does the 'usual' normalisation,           */
  /* the factor  2.0                  accounts for 1-sidedness of the spectrum, */
  /* the factor  1/nsegment           does the averaging over segments.         */
  for (i=0; i<DF->FTSize; ++i) {
    DF->powspec[i] = log(DF->powspec[i]) + psdcoef;
  }
  /*     /!\   'winss' coefficient not (yet?) considered !!    */
  if (verbose) printf(" | %d data segments averaged for PSD estimation.\n", nsegment);
  DF->PsdEstimateN = nsegment;
  /* up to now, "DF->PsdEstimateN" denoted the (maximum) desired number of data   */
  /* segments to average over, now it gives the number of actually used segments. */
  return PsdOk;
}


void simulatePsdEstimation(DataFramework *DF)
/* simulates the process of PSD-estimation for "real" data, but               */
/* using simulated noise... may in fact be more accurate than just            */
/* plugging in the true spectrum, as this will also account for the effects   */
/* of windowing, and anyway more closely resemble the algorithm's performance */
/* for actual data.                                                           */
{
  int i;
  double *noise;
  double complex *noiseFT;
  long j;
  double power, psdcoef;

  if (verbose) printf(" | simulating noise PSD estimation (%s)\n", DF->ifo->name);

  noise   = (double*) malloc(sizeof(double)*DF->dataSize);
  noiseFT = (double complex*) malloc(sizeof(double complex) * (DF->dataSize/2 + 1));
  for (j=0; j<DF->FTSize; ++j)
    DF->powspec[j] = 0.0;

  for (i=0; i<DF->PsdEstimateN; ++i){
    generateNoise(DF->dataDeltaT, DF->dataSize, noise, DF->noiseModel);
    FTexec(DF, noise, noiseFT);
    for (j=0; j<DF->FTSize; ++j){
      power = cabs(noiseFT[j]);
      power *= power;
      DF->powspec[j] += power;
    }
  } 
  free(noise);
  free(noiseFT);
  psdcoef = log(2.0) + log(DF->dataDeltaT) - log(DF->dataSize) - log(DF->PsdEstimateN);
  /* the factor  dataDeltaT/dataSize  does the 'usual' normalisation,           */
  /* the factor  2.0                  accounts for 1-sidedness of the spectrum, */
  /* the factor  1/DF->PsdEstimateN   does the averaging over segments.         */
  for (i=0; i<DF->FTSize; ++i) {
    DF->powspec[i] = log(DF->powspec[i]) + psdcoef;
  }
}


double *downsample(double data[], int *datalength, int factor)
/* `datalength' returns the length of the new data vector                                 */
{
  static LALStatus status;
  REAL4TimeSeries  lalTS;
  ResampleTSParams resampPars;
  int i, newlength;
  double *decimated;

  /* little sanity check: */
  if ((factor!=1) && (factor!=2) && (factor!=4) && (factor!=8))
    printf(" : ERROR: only downsampling factors 1, 2, 4, or 8 allowed!\n");
  if (factor > 1) {
    /* first downsample time series                                      */
    /* (this part is basically copied from "ResampleTimeSeriesTest.c"):  */
    memset(&lalTS, 0, sizeof(REAL4TimeSeries));
    LALSCreateVector(&status, &(lalTS.data), *datalength);
    lalTS.sampleUnits = lalADCCountUnit;
    lalTS.deltaT = 1.0 / ((double) factor);
    lalTS.epoch.gpsSeconds = 0;
    resampPars.deltaT = 1.0;
    resampPars.filterType = LDASfirLP;
    /*  (the `LDASfirLP' option is used here because it  */ 
    /*  guarantees a fixed number of corrupted samples)  */
    for (i=0; i<lalTS.data->length; ++i)
      lalTS.data->data[i] = data[i];
    LALResampleREAL4TimeSeries(&status, &lalTS, &resampPars);
    /* according to LAL documentation,                                  */
    /* there should now be  20  corrupted (=zero) samples at both ends: */
    newlength = (*datalength - 2*20*factor) / factor;
    if (newlength != (lalTS.data->length - 40))
      printf(" : ERROR: vector size mismatch in `downsample()': newlength=%d vs. data->length-40=%d\n",
	     newlength, lalTS.data->length - 40);
    decimated = (double*) malloc(newlength * sizeof(double));
    for (i=0; i<(lalTS.data->length - 40); ++i)
      decimated[i] = lalTS.data->data[i+20];
    LALSDestroyVector(&status, &lalTS.data);
  }
  else {
    newlength = *datalength;
    decimated = (double*) malloc(newlength * sizeof(double));
    for (i=0; i<(*datalength); ++i)
      decimated[i] = data[i];
  }

  *datalength = newlength;
  return(decimated);
}


double tukeywindow(int j, int N, double r)
/* Tukey window... for r=0 equal to rectangular window, for r=1 equal to Hann window. */
/* 0 < r < 1 denotes the fraction of the window in which it behaves sinusoidal.       */
/* j = 0, 1, ..., N-1.                                                                */
{
  double win = 1.0;
  if (j>(((double)N)/2.0)) j = N-j;
  if (((double)j) < (r*(((double)N)/2.0)))
    win = 0.5*(1.0-cos(((2.0*pi)/r)*(((double)j)/((double)N))));
  return win;
}


void FTexec(DataFramework *DF, double *input, fftw_complex *output)
/*******************************************/
/*  Executes a Fourier transform based on  */
/*  the `initDF()' settings given in `DF'. */
/*  `input'  is assumed to be a vector     */
/*           of length `DF->dataSize'.     */
/*  `output' is assumed to be a            */
/*           pre-allocated (!!) vector     */
/*           of length `fw->FTN'.          */
/*  `input' is NOT overwritten during FT.  */
/* * * * * * * * * * * * * * * * * * * * * */
/*            Not too efficient            */
/*      (due to back-and-forth-copying     */
/*           of input and output)          */
/*            but simple instead.          */
/*******************************************/
{
  long i;
  if (output==NULL) 
    printf("  !! ERROR: unallocated output vector in `FTexec()' !!\n");
  /* copy AND WINDOW (!) the input data: */
  for (i=0; i<DF->dataSize; ++i)
    DF->FTin[i] = input[i] * DF->window[i];
  fftw_execute(DF->FTplan);
  for (i=0; i<DF->FTSize; ++i)
    output[i] = DF->FTout[i];
}


double spectrumIniLigo(double f)
/* (logarithmic) one-sided power spectral density of initial-LIGO noise */
/* following Damour/Iyer/Sathyaprakash (2001), PRD 63 044023; table iv  */
/* and LAL's "LALLIGOIPsd()" function.                                  */
/* NOTE: in this implementation, the PSD is limited above               */
/* at 2e10 times its "sweet spot" value at 150 Hz.                      */
{
  double logpsd150, upperlimit, logpsd, x;
  if (f>0.0) {
    logpsd150 = log(9e-46) + log(pow(4.49,-56.0) + 0.16 + 0.52 + 0.32);
    upperlimit = logpsd150 + log(2e10);
    x = f/150.0;
    logpsd = log(9e-46) + log(pow(4.49*x,-56.0) + 0.16*pow(x,-4.52) + 0.52 + 0.32*x*x);
    if (logpsd>upperlimit)
      logpsd = upperlimit;
  }
  else logpsd = -HUGE_VAL;
  return logpsd;
}


double spectrumAdvLigo(double f)
/* (logarithmic) one-sided power spectral density of advanced-LIGO noise */
/* following LAL's "LALAdvLIGOPsd()" function.                           */
/* NOTE: in this implementation, the PSD is limited above                */
/* at 2e10 times its "sweet spot" value at 215 Hz.                       */
{
  double logpsd215, upperlimit, logpsd, x, x2;
  if (f>0.0) {
    logpsd215 = log(1e-49) + log(1.0 - 5.0 + 111.0*0.5/1.5);
    upperlimit = logpsd215 + log(2e10);
    x = f/215;
    x2 = x*x;
    logpsd = log(1e-49) + log(pow(x,-4.14) - 5.0/x2 + 111.0*(1.0-x2+0.5*x2*x2)/(1.0+0.5*x2));
    if (logpsd>upperlimit)
      logpsd = upperlimit;
  }
  else logpsd = -HUGE_VAL;
  return logpsd;
}


double spectrumVirgo(double f)
/* (logarithmic) one-sided power spectral density of (initial?) Virgo noise */
/* following LAL's "LALVIRGOPsd()" function.                                */
/* NOTE: in this implementation, the PSD is limited above                   */
/* at 2e10 times its "sweet spot" value at 281 Hz.                          */
{
  double logpsd281, upperlimit, logpsd, x;
  if (f>0.0) {
    x = 281.0/500.0;
    logpsd281 = log(10.2e-46) + log(pow(7.87*x,-4.8) + 6.0/17.0/x + 1.0 + x*x);
    upperlimit = logpsd281 + log(2e10);
    x = f/500.0;
    logpsd = log(10.2e-46) + log(pow(7.87*x,-4.8) + 6.0/17.0/x + 1.0 + x*x);
    if (logpsd>upperlimit)
      logpsd = upperlimit;
  }
  else logpsd = -HUGE_VAL;
  return logpsd;
}


double spectrumGeo(double f)
/* (logarithmic) one-sided power spectral density of GEO noise          */
/* following Damour/Iyer/Sathyaprakash (2001), PRD 63 044023; table iv  */
/* NOTE: in this implementation, the PSD is limited above               */
/* at 2e10 times its "sweet spot" value at 150 Hz.                      */
{
  double logpsd150, upperlimit, logpsd, x, x2;
  if (f>0.0) {
    logpsd150 = -46.0*log(10.0) + log(pow(3.4,-30.0) + 34.0 + 10.0/1.5);
    upperlimit = logpsd150 + log(2e10);
    x  = f/150.0;
    x2 = x*x;
    logpsd = -46.0*log(10.0) + log(pow(3.4*x,-30.0) + 34.0/x + 20.0*(1-x2+0.5*x2*x2)/(1+x2/2.0));
    if (logpsd>upperlimit)
      logpsd = upperlimit;
  }
  else logpsd = -HUGE_VAL;
  return logpsd;
}


void normalise(double vec[3])
/*   vec  :=  vec / |vec|   */
{
  double length = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
  vec[0] /= length;
  vec[1] /= length;
  vec[2] /= length;
}


void rotate(double x[3], double angle, double axis[3])
/* rotates vector x clockwise around `axis'               */
/* (looking along axis while it is pointing towards you). */
/*   !!  `axis' must be a UNIT VECTOR  !!                 */
{
  int i, j;
  double cosa = cos(-angle);
  double sina = sin(-angle);
  double R[3][3] = {{cosa+axis[0]*axis[0]*(1.0-cosa), 
                     axis[0]*axis[1]*(1.0-cosa)-axis[2]*sina,
                     axis[0]*axis[2]*(1.0-cosa)+axis[1]*sina},
                    {axis[1]*axis[0]*(1.0-cosa)+axis[2]*sina,
                     cosa+axis[1]*axis[1]*(1.0-cosa),
                     axis[1]*axis[2]*(1.0-cosa)-axis[0]*sina},
                    {axis[2]*axis[0]*(1.0-cosa)-axis[1]*sina,
                     axis[2]*axis[1]*(1.0-cosa)+axis[0]*sina,
                     cosa+axis[2]*axis[2]*(1.0-cosa)}};
  double result[3] = {0.0, 0.0, 0.0};
  for (i=0; i<3; ++i)
    for (j=0; j<3; ++j)
      result[i] += R[i][j]*x[j];
  for (i=0; i<3; ++i) x[i] = result[i];
}


int righthanded(double x[3], double y[3], double z[3])
/* Determines whether vectors x,y & z constitute a right-handed system */
/* by checking the sign of the triple product or det(x,y,z).           */
{
  double spatprodukt =   x[0]*y[1]*z[2] + y[0]*z[1]*x[2] + z[0]*x[1]*y[2]
                       - z[0]*y[1]*x[2] - x[0]*z[1]*y[2] - y[0]*x[1]*z[2];
  int result = (spatprodukt >= 0.0) ? 1 : 0;
  return result;
}


void orthoproject(double x[3], double vec1[3], double vec2[3])
/* Determines the orthogonal projection of vector x onto the span of */
/* the two ORTHONORMAL (!) vectors vec1 and vec2.                    */
{
  double sprod1 = x[0]*vec1[0] + x[1]*vec1[1] + x[2]*vec1[2];
  double sprod2 = x[0]*vec2[0] + x[1]*vec2[1] + x[2]*vec2[2];
  x[0] = sprod1*vec1[0] + sprod2*vec2[0];
  x[1] = sprod1*vec1[1] + sprod2*vec2[1];
  x[2] = sprod1*vec1[2] + sprod2*vec2[2];
}


double angle(double x[3], double y[3])
/* Determines the angle between vectors x & y. */
{
  double sprod = x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
  double absx  = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  double absy  = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
  return acos(sprod/(absx*absy));
}


void coord2vec(double lati, double longi, double x[3])
/* Turns geographical coordinates (latitude & longitude) into a vector -     */
/* Result is a unit (!) vector referring to the (right-handed) coordinate    */
/* system spanned by the three vectors pointing from geocenter to:           */
/*   x) intersection of greenwich meridian and equator                       */
/*   y) intersection of 90 deg. East meridian and equator                    */
/*   z) north pole                                                           */
{
  double sinlati = sin(lati);
  double coslati = sqrt(1.0-sinlati*sinlati);
  x[0] = cos(longi) * coslati;  /* `Greenwich' dimension  */
  x[1] = sin(longi) * coslati;  /* `Ganges' dimension     */
  x[2] = sinlati;               /* `North Pole' dimension */
}


void vec2coord(double x[3], double *lati, double *longi)
/* Derive geographical coordinates from a vector (see also previous function). */
{
  double greenwich[3] = {1.0, 0.0, 0.0};
  double ganges[3]    = {0.0, 1.0, 0.0};
  double northpole[3] = {0.0, 0.0, 1.0};
  double dummy[3]     = {x[0], x[1], x[2]};
  *lati = 0.5*pi - angle(northpole, x);
  orthoproject(dummy, greenwich, ganges);
  *longi = angle(greenwich, dummy);
  if (righthanded(greenwich,northpole,dummy))
    *longi *= -1.0;
}


void ifoInitComputations(interferometer *ifolist, int ifoN)
/* Compute location/orientation vectors            */
/* associated with interferometers and             */
/* derived from given location/orientation angles. */ 
/* Assumes angles to be specified already          */
/* and is called within "ifoInit()", see below.    */
{
  int i;
  double merinormal[3];
  double eccentricity2, curvatureradius;
  for (i=0; i<ifoN; ++i){
    /* place arms on equator plane, so that its designated N-S-direction */
    /* is aligned with its meridian plane:                               */
    ifolist[i].rightArmVector[0] = -cos(ifolist[i].longitude + ifolist[i].rightArmAngle);
    ifolist[i].rightArmVector[1] = -sin(ifolist[i].longitude + ifolist[i].rightArmAngle);
    ifolist[i].rightArmVector[2] = 0.0;
    ifolist[i].leftArmVector[0]  = -cos(ifolist[i].longitude + ifolist[i].rightArmAngle + 0.5*pi);
    ifolist[i].leftArmVector[1]  = -sin(ifolist[i].longitude + ifolist[i].rightArmAngle + 0.5*pi);
    ifolist[i].leftArmVector[2]  = 0.0;
    ifolist[i].normalVector[0]   = 0.0;
    ifolist[i].normalVector[1]   = 0.0;
    ifolist[i].normalVector[2]   = 1.0;
    /* Determine normal vector of meridian plane: */
    merinormal[0] = cos(ifolist[i].longitude - 0.5*pi);
    merinormal[1] = sin(ifolist[i].longitude - 0.5*pi);
    merinormal[2] = 0.0;
    /* The three vectors:                                                      */
    /*   x) from geocenter to intersection of ifo meridian with equator plane  */
    /*   y) from geocenter to north pole                                       */
    /*   z) the above normal vector                                            */
    /* again form another (orthonormal) right-handed system.                   */
    /* Now turn all arms clockwise around the normal vector of meridian plane: */
    rotate(ifolist[i].rightArmVector, pi/2.0 - ifolist[i].latitude, merinormal);
    rotate(ifolist[i].leftArmVector,  pi/2.0 - ifolist[i].latitude, merinormal);
    rotate(ifolist[i].normalVector,   pi/2.0 - ifolist[i].latitude, merinormal);
    /* initialise the ifo position (!NOT! unit-) vector: */
    coord2vec(ifolist[i].latitude, ifolist[i].longitude, ifolist[i].positionVector);
    eccentricity2  = earthFlattening*(2.0-earthFlattening);  /* (squared eccentricity) */
    curvatureradius = earthRadiusEquator/sqrt(1.0-eccentricity2*pow(sin(ifolist[i].latitude),2.0));
    ifolist[i].positionVector[0] *= curvatureradius;
    ifolist[i].positionVector[1] *= curvatureradius;
    ifolist[i].positionVector[2] *= curvatureradius*(1.0-eccentricity2);
  }
}


void ifoInit(interferometer *ifolist[], int *ifoN)
/* this function fills up the "database" of interferometers */
/* (names, locations, orientations)                         */
/* and also derives vectors etc.                            */
/* Data is taken from:                                      */
/*   Gravitational wave detector sites                      */
/*   http://arxiv.org/abs/gr-qc/9607075                     */
{
  int i;
  interferometer *ilist=NULL;
  *ifoN = 4;  /* for now 4 Ifos: Hanford, Livingston, Pisa, Hannover. */
  ilist = (interferometer*) malloc(sizeof(interferometer) * (*ifoN));
  for (i=0; i<(*ifoN); ++i)
    ilist[i].name = (char*) malloc(sizeof(char)*32); /* for now limited to 31 characters */
  /* fill in HANFORD details:    */
  strcpy(ilist[0].name, "LIGO-Hanford");
  ilist[0].latitude      = (  46.45/180.0)*pi;
  ilist[0].longitude     = (-119.41/180.0)*pi;
  ilist[0].rightArmAngle = (  36.80/180.0)*pi;
  /* fill in LIVINGSTON details: */
  strcpy(ilist[1].name, "LIGO-Livingston");
  ilist[1].latitude      = (  30.56/180.0)*pi;
  ilist[1].longitude     = ( -90.77/180.0)*pi;
  ilist[1].rightArmAngle = ( 108.00/180.0)*pi;
  /* fill in PISA details:       */
  strcpy(ilist[2].name, "VIRGO-Pisa");
  ilist[2].latitude      = (  43.63/180.0)*pi;
  ilist[2].longitude     = (  10.50/180.0)*pi;
  ilist[2].rightArmAngle = ( 341.50/180.0)*pi;
  /* fill in HANNOVER details:   */
  strcpy(ilist[3].name, "GEO-Hannover");
  ilist[3].latitude      = (  52.25/180.0)*pi;
  ilist[3].longitude     = (   9.81/180.0)*pi;
  ilist[3].rightArmAngle = ((291.61+2.165)/180.0)*pi;

  /* now compute vectors etc: */
  ifoInitComputations(ilist, *ifoN);

  *ifolist = ilist;
}


void clearIfo(interferometer *ifolist, int *ifoN)
{
  int i;
  if (ifolist != NULL) {
    for (i=0; i<*ifoN; ++i) free(ifolist[i].name);
    free(ifolist);
  }
  ifolist = NULL;
  *ifoN = 0;
}


void localParameters(vector *parameter, interferometer *ifo,
                     double *timeshift, double *polarisation, double *altitude, double *azimuth)
/* determine the "local" parameters, with respect to the particular interferometer `ifo',  */
/* and derived from interferometer's location/orientation and line-of-sight.               */
{
  double lineofsight[3];
  double dummyvec[3];
  double polarvec[3];
  double scalprod;
  int i;
  coord2vec(vectorGetValue(parameter,"latitude"), vectorGetValue(parameter,"longitude"),
	    lineofsight);
  /*-- determine COALESCENCE TIME  (with respect to ifo'): --*/
  scalprod =   ifo->positionVector[0] * lineofsight[0]
             + ifo->positionVector[1] * lineofsight[1]
             + ifo->positionVector[2] * lineofsight[2];
  *timeshift = (-1.0) * scalprod / c;
  /* (negative timeshift means signal arrives earlier than at geocenter etc.) */
  /*-- determine ALTITUDE (with respect to ifo'): --*/
  *altitude = angle(ifo->normalVector, lineofsight); 
  /*-- determine AZIMUTH (w.r.t. ifo'):           --*/
  /* project line of sight into ifo' arm plane: */
  for (i=0; i<3; ++i) dummyvec[i] = lineofsight[i];
  orthoproject(dummyvec, ifo->rightArmVector, ifo->leftArmVector);
  *azimuth = angle(dummyvec, ifo->rightArmVector);
  if (!righthanded(ifo->rightArmVector, dummyvec, ifo->normalVector))
    *azimuth = 2.0*pi - *azimuth;
  /*-- determine POLARISATION (w.r.t. ifo'):      --*/
  /* construct a `polarisation vector': */
  scalprod = lineofsight[2]; /* product of north pole vector (0,0,1) and `lineofsight' */
  polarvec[0] = 0.0 - scalprod*lineofsight[0];
  polarvec[1] = 0.0 - scalprod*lineofsight[1];
  polarvec[2] = 1.0 - scalprod*lineofsight[2];
  /* `polarvec' by now is the north vector tilted so it is orthogonal to line-of-sight */
  rotate(polarvec, vectorGetValue(parameter,"polarisation") - 0.5*pi, lineofsight);
  /* `polarvec' now is the normal vector of the polarisation plane                 */
  /* (the plane spanned by line-of-sight and polarisation).                        */
  /* Line-of-sight, polarisation and `polarvec' form a right-handed system.        */
  /* Now rotate `dummyvec' so it is orthogonal to line-of-sight/normalvec plane:   */
  rotate(dummyvec, 0.5*pi, ifo->normalVector);
  /* Line-of-sight, `normalvec' and `dummyvec' form a right-handed system as well. */
  *polarisation = angle(dummyvec, polarvec);
  if (righthanded(lineofsight, dummyvec, polarvec))
    *polarisation = pi - *polarisation;
  /* -so polarisation measures counterclockwise angle looking along `lineofsight'. */
}


void antennaepattern(double altitude, double azimuth, double polarisation,
                     double *Fplus, double *Fcross)
/*  Antennae pattern (or beam-pattern) functions,  */
/*  following Blanchet (2001) gr-qc/0104084        */
{
  double cosalti   = cos(altitude);
  double sin2azi   = sin(2.0*azimuth);
  double cos2azi   = cos(2.0*azimuth);
  double sin2polar = sin(2.0*polarisation);
  double cos2polar = cos(2.0*polarisation);
  *Fplus  = 0.5*(1.0+cosalti*cosalti)*cos2azi*cos2polar
            + cosalti*sin2azi*sin2polar;                  /* (3.11) */
  *Fcross = -0.5*(1.0+cosalti*cosalti)*cos2azi*sin2polar
            + cosalti*sin2azi*cos2polar;                  /* (3.12) */
}


double mc2mass1(double mc, double eta)
/* mass 1 (the smaller one) for given mass ratio & chirp mass */
{
  double root = sqrt(0.25-eta);
  double fraction = (0.5+root) / (0.5-root);
  return mc * (pow(1+fraction,0.2) / pow(fraction,0.6));
}


double mc2mass2(double mc, double eta)
/* mass 2 (the greater one) for given mass ratio & chirp mass */
{
  double root = sqrt(0.25-eta);
  double inversefraction = (0.5-root) / (0.5+root);
  return mc * (pow(1+inversefraction,0.2) / pow(inversefraction,0.6));
}


double mc2mt(double mc, double eta)
/* total mass (mt) for given mass ratio & chirp mass */
{
  double root = sqrt(0.25-eta);
  double fraction = (0.5+root) / (0.5-root);
  double inversefraction = (0.5-root) / (0.5+root);
  return mc * ((pow(1+fraction,0.2) / pow(fraction,0.6)) 
               + (pow(1+inversefraction,0.2) / pow(inversefraction,0.6)));
}


double logJacobianMcEta(double mc, double eta)
/* posterior multiplier for transformed parameters */
/* (jacobian of inverse tranformation)             */
/* (mc & eta  instead of  m1 & m2)                 */
{
  double result;
  double m1, m2, msum, mprod;
  double term1, term2, term3, term4;
  double a,b,c,d;
  /* Factor for the  (m1,m2) --> (mc,eta)  transform: */
  msum  = mc2mt(mc,eta) * Msun;
  m1    = mc2mass1(mc, eta) * Msun;
  m2    = msum-m1;
  mprod = m1*m2;
  term1 = 0.6*pow(msum,-0.2);
  term2 = 0.2*pow(msum,-1.2)*pow(mprod,0.6);
  term3 = pow(msum,-2.0);
  term4 = 2*mprod*pow(msum,-3.0);
  a = pow(m1,-0.4)*pow(m2,0.6)*term1 - term2;
  b = m2*term3-term4;
  c = pow(m1,0.6)*pow(m2,-0.4)*term1 - term2;
  d = m1*term3-term4;
  result =  -log(fabs(a*d - b*c));
  return result;
}


void signaltemplate(DataFramework *DF, int waveform, vector *parameter, double complex *output)
/* 'wrapper' function for signal waveform templates.                                */
/* Checks the "waveform" character string and calls                                 */
/* the corresponding waveform-generating function (with "parameter" argument)       */
/* and returns the frequancy-domain template in the (pre-allocated, double complex) */
/* "output" vector of length "DF->FTSize".                                          */
{
  double Fplus, Fcross;
  antennaepattern(vectorGetValue(parameter,"altitude"), 
                  vectorGetValue(parameter,"azimuth"), 
                  vectorGetValue(parameter,"polarisation"),
                  &Fplus, &Fcross);
  if (waveform == iR2PN)               /* restricted 2.0 PN              */
    templateR2PN(DF, parameter, Fplus, Fcross, output);
  else if (waveform == i20SP)          /* 2.0 PN stationary phase        */
    templateStatPhase(DF, parameter, Fplus, Fcross, output, 2.0);
  else if (waveform == i25SP)          /* 2.5 PN stationary phase        */
    templateStatPhase(DF, parameter, Fplus, Fcross, output, 2.5);
  else if (waveform == i2025)          /* 2.0 PN amplitude, 2.5 PN phase */
    template2025(DF, parameter, Fplus, Fcross, output);
  else if (waveform == i2535)          /* 2.5 PN amplitude, 3.5 PN phase */
    template2535(DF, parameter, Fplus, Fcross, output);
  /*-- LAL templates... --*/
  else if (waveform == iLALTT2PN00)    /* LAL Taylor T2 Newtonian        */
    templateLAL(DF, parameter, Fplus, Fcross, output, TaylorT2, newtonian);
  else if (waveform == iLALTT2PN10)    /* LAL Taylor T2 1.0PN            */
    templateLAL(DF, parameter, Fplus, Fcross, output, TaylorT2, onePN);
  else if (waveform == iLALTT2PN15)    /* LAL Taylor T2 1.5PN            */
    templateLAL(DF, parameter, Fplus, Fcross, output, TaylorT2, onePointFivePN);
  else if (waveform == iLALTT2PN20)    /* LAL Taylor T2 2.0PN            */
    templateLAL(DF, parameter, Fplus, Fcross, output, TaylorT2, twoPN);
  else if (waveform == iLALTT3PN00)    /* LAL Taylor T3 Newtonian        */
    templateLAL(DF, parameter, Fplus, Fcross, output, TaylorT3, newtonian);
  else if (waveform == iLALTT3PN10)    /* LAL Taylor T3 1.0PN            */
    templateLAL(DF, parameter, Fplus, Fcross, output, TaylorT3, onePN);
  else if (waveform == iLALTT3PN15)    /* LAL Taylor T3 1.5PN            */
    templateLAL(DF, parameter, Fplus, Fcross, output, TaylorT3, onePointFivePN);
  else if (waveform == iLALTT3PN20)    /* LAL Taylor T3 2.0PN            */
    templateLAL(DF, parameter, Fplus, Fcross, output, TaylorT3, twoPN);
  else if (waveform == iLALIMRPhenomA) /* LAL Phenomenological           */
    templateLAL(DF, parameter, Fplus, Fcross, output, IMRPhenomA, pseudoFourPN);
  else if (waveform == iLALEOBNR)      /* LAL EOBNR                      */
    templateLAL(DF, parameter, Fplus, Fcross, output, EOBNR, pseudoFourPN);
  /* burst: */
  else if (waveform == bSineGaussian)  /* Sine-Gaussian burst            */
    templateSineGaussianBurst(DF, parameter, Fplus, Fcross, output);
  else printf(" : ERROR: requesting undefined template waveform (%d) in 'signaltemplate()'!\n", 
              waveform);
}


void template2025(DataFramework *DF, vector *parameter, double Fplus, double Fcross,
                  double complex *output)
/*************************************************************/
/* returns the (numerically FT'd) frequency-domain template. */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Time domain template computation following                */
/* Blanchet (2001), gr-qc/0104084.                           */
/* The approximation is:   2.0 PN in amplitude  and          */
/*                         2.5 PN in phase.                  */
/* (formula numbers (x.xx) refer to the Blanchet paper)      */
/*************************************************************/
{
  double eta      = vectorGetValue(parameter,"massratio");
  double m1       = mc2mass1(vectorGetValue(parameter,"chirpmass"), eta); /* (in units of Ms) */
  double m2       = mc2mass2(vectorGetValue(parameter,"chirpmass"), eta);
  double mt       = m1 + m2;
  double dmm      = (m2 - m1) / mt; /*  = (delta m) / mt  (dimensionless) */
  double log_mt   = log(mt) + log(Msun);  /* (in Kg) */
  double log_eta  = log(eta);
  double eta2     = eta * eta;
  double log_mu   = log_eta + log_mt;
  double t, phi, psi, omega, oldomega=0.0;
  double phicoef[5]   = {-1.0/eta,                                                /*  (6.13)  */
                         0.460689484126984 + 0.5729166666666666*eta,
                         -2.356194490192345,
                         0.6418722070533943 + (1.103961278521825 + 0.90576171875*eta)*eta,
                         -0.7057224708076262 - 0.02300971181828462*eta};
  double omegacoef[5] = {exp(79.90678935048672 - log_mt),                         /*  (6.14)  */
                         0.2764136904761905 + 0.34375 * eta,
                         -0.9424777960769379,
                         0.1283744414106789 + (0.2207922557043651 + 0.18115234375*eta)*eta,
                         -1.129155953292202 - 0.03681553890925539*eta};
  double psicoef = exp(-81.29308371160661 + log_mt);                              /*  (6.12)  */
  double taucoef = 80.37679297973246 + log_eta - log_mt;                          /*  (4.17)  */
  double log_tau, tau18, tau28, tau38, tau58, tau68, tau78, tau88;
  double ci  = cos(vectorGetValue(parameter,"inclination"));
  double ci2 = ci*ci,     ci4 = ci2*ci2,   ci6 = ci4*ci2;
  double si2 = (1.0-ci2), si  = sqrt(si2), si4 = si2*si2;
  double h_plus, h_cross;
  double Hp00, Hp05, Hp10, Hp15, Hp20;
  double Hc00, Hc05, Hc10, Hc15, Hc20;
  double sin1psi, sin2psi, sin3psi, sin4psi, sin5psi, sin6psi;
  double cos1psi, cos2psi, cos3psi, cos4psi, cos5psi, cos6psi;
  double constfactor=exp(-61.77448272506146+log_mu-(vectorGetValue(parameter,"logdistance")+log(Mpc))); /*(6.01)*/
  double sqrtx, xcoef = -81.98623089216656+log_mt;   /*  = log((G*mt)/(c^3))          (6.01)  */
     int i;
  double *timedomainwaveform=NULL;
  double tc = vectorGetValue(parameter,"time");
  double phase = vectorGetValue(parameter,"phase");

  /* fill `timedomainwaveform' with time-domain template: */
  timedomainwaveform = (double*) malloc(sizeof(double)*DF->dataSize);
  for (i=0; i<DF->dataSize; ++i){
    t = (tc - DF->dataStart) - ((double)i)*DF->dataDeltaT;  /* (time to t_c) = "(t_c-t)" in (4.17) */
    /* == t_c - (tstart + i*deltaT) */
    if (t>0.0) {
      log_tau = taucoef + log(t);                                 /*                  (4.17)  */
      tau28   = exp(0.25 * log_tau);                              /* = tau ^ (2/8)            */
      tau38   = exp(0.375 * log_tau);                             /* = tau ^ (3/8)            */
      tau58   = tau28 * tau38;                                    /* = tau ^ (5/8)            */
      tau68   = tau38 * tau38;                                    /* = tau ^ (6/8)            */
      tau78   = tau28 * tau58;                                    /* = tau ^ (7/8)            */
      tau88   = tau28 * tau68;                                    /* = tau ^ (8/8)            */
      tau18   = tau38 / tau28;                                    /* = tau ^ (1/8)            */
      omega   = omegacoef[0] * (1.0/tau38 + omegacoef[1]/tau58 + omegacoef[2]/tau68
				+ omegacoef[3]/tau78 + omegacoef[4]/tau88);       /*  (6.14)  */
      if (!(omega >= oldomega)){   /*  (frequency starts decreasing)  */
        timedomainwaveform[i] = 0.0; 
        omega = oldomega;
      }
      else {
        oldomega= omega;
        sqrtx   = exp(0.5 * (2.0/3.0) * (xcoef+log(omega)));        /* = x ^ (1/2)    (4.13)  */
	phi     = phase + phicoef[0]*(tau58 + phicoef[1]*tau38 + phicoef[2]*tau28
                        + phicoef[3]*tau18 + phicoef[4]*log_tau);                 /*  (6.13)  */
	psi     = phi - psicoef * omega * (log(omega)-2.754167798283500);         /*  (6.12)  */
	sin1psi = sin(psi);      cos1psi = cos(psi);
	sin2psi = sin(2.0*psi);  cos2psi = cos(2.0*psi);
	sin3psi = sin(3.0*psi);  cos3psi = cos(3.0*psi);
	sin4psi = sin(4.0*psi);  cos4psi = cos(4.0*psi);
	sin5psi = sin(5.0*psi);  cos5psi = cos(5.0*psi);
	sin6psi = sin(6.0*psi);  cos6psi = cos(6.0*psi);
	Hp00    = -(1.0+ci2)*cos2psi;                                             /*  (6.02)  */
	Hp05    = -(si/8.0)*dmm * ((5.0+ci2)*cos1psi - 9.0*(1.0+ci2)*cos3psi);    /*  (6.03)  */
	Hp10    = 0.1666666666666667*((19.0+9.0*ci2-2.0*ci4)-eta*(19.0-11.0*ci2-6.0*ci4))*cos2psi
	          - 1.333333333333333*si2*(1.0+ci2)*(1.0-3.0*eta)*cos4psi;        /*  (6.04)  */
	Hp15    = (si/192.0)*dmm * (((57.0 + 60.0*ci2-ci4) - 2.0*eta*(49.0-12.0*ci2-ci4))*cos1psi
	          -13.5*((73.0+40.0*ci2-9.0*ci4) - 2.0*eta*(25.0-8.0*ci2-9.0*ci4))*cos3psi
	          +312.5*(1.0-2.0*eta)*si2*(1.0+ci2)*cos5psi) - 6.283185307179586*(1.0+ci2)*cos2psi;
	                                                                          /*  (6.05)  */
	Hp20    = 0.008333333333333333*((22.0+396.0*ci2+145.0*ci4-5.0*ci6)        
	          + 1.666666666666667*eta*(706.0-216.0*ci2-251.0*ci4+15.0*ci6)
	          -5.0*eta2*(98.0-108.0*ci2+7.0*ci4+5.0*ci6))*cos2psi
	          +0.1333333333333333*si2*((59.0+35.0*ci2-8.0*ci4)
	          -1.666666666666667*eta*(131.0+59.0*ci2-24.0*ci4)
	          +5.0*eta2*(21.0-3.0*ci2-8.0*ci4))*cos4psi
	          -2.025*(1.0-5.0*eta+5.0*eta2)*si4*(1.0+ci2)*cos6psi
	          +si/40.0*dmm*((11.0+7.0*ci2+10.0*(5.0+ci2)*0.6931471805599453)*sin1psi 
	          -15.70796326794897*(5.0+ci2)*cos1psi
	          -79.52442081079560*(1.0+ci2)*sin3psi+424.115008234622*(1.0+ci2)*cos3psi);
	                                                                          /*  (6.06)  */
	Hc00    = -2.0*ci*sin2psi;                                                /*  (6.07)  */
	Hc05    = -0.75*si*ci*dmm*(sin1psi-3.0*sin3psi);                          /*  (6.08)  */
	Hc10    = (ci/3.0)*((17.0-4.0*ci2)-eta*(13.0-12.0*ci2))*sin2psi
	          -2.666666666666667*(1.0-3.0*eta)*ci*si2*sin4psi;                /*  (6.09)  */
	Hc15    = ((si*ci)/96.0)*dmm*(((63.0-5.0*ci2)-2.0*eta*(23.0-5.0*ci2))*sin1psi
	          -13.5*((67.0-15.0*ci2)-2.0*eta*(19.0-15.0*ci2))*sin3psi
	          +312.5*(1.0-2.0*eta)*si2*sin5psi)-12.56637061435917*ci*sin2psi; /*  (6.10)  */
	Hc20    = (ci/60.0)*((68.0+226.0*ci2-15.0*ci4)+1.666666666666667*eta*(572.0-490.0*ci2+45.0*ci4)
	          -5.0*eta2*(56.0-70.0*ci2+15.0*ci4))*sin2psi
	          +0.2666666666666667*ci*si2*((55.0-12.0*ci2)-1.666666666666667*eta*(119.0-36.0*ci2)
	          +5.0*eta2*(17.0-12.0*ci2))*sin4psi
	          -4.05*(1.0-5.0*eta+5.0*eta2)*ci*si4*sin6psi
	          -0.15*si*ci*dmm*(9.931471805599454*cos1psi+15.70796326794897*sin1psi
	          -26.50814027026520*cos3psi - 141.3716694115407*sin3psi);        /*  (6.11)  */
	h_plus  = h_cross = constfactor * sqrtx*sqrtx;
	h_plus  *= Hp00 + sqrtx*(Hp05 + sqrtx*(Hp10 + sqrtx*(Hp15 + sqrtx*Hp20)));
	h_cross *= Hc00 + sqrtx*(Hc05 + sqrtx*(Hc10 + sqrtx*(Hc15 + sqrtx*Hc20)));/*  (6.01)  */
	timedomainwaveform[i] = Fplus*h_plus + Fcross*h_cross;                    /*  (3.10)  */
      }
    }
    else timedomainwaveform[i] = 0.0; 
  }
  /* FT'd waveform is returned in "output": */
  FTexec(DF, timedomainwaveform, output);
  free(timedomainwaveform);
}


void template2535(DataFramework *DF, vector *parameter, double Fplus, double Fcross,
                  double complex *output)
/***************************************************************************/
/* returns the (numerically FT'd) frequency-domain template.               */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* 2.5PN amplitude / 3.5PN phase  time domain template computation         */
/* following Blanchet (2001), gr-qc/0105099                                */
/* and Arun et al. (2004), CQG 21(15):3771-3801                            */
/* Formula numbers (x.xx) refer to the `old' Blanchet paper gr-qc/0104084, */
/* numbers (xx) refer to the more recent gr-qc/0105099.                    */
/* Numbers referring to Arun et al (2004) are explicitly marked.           */
/***************************************************************************/
{
  double eta        = vectorGetValue(parameter,"massratio");
  double m1         = mc2mass1(vectorGetValue(parameter,"chirpmass"), eta); /* (in units of Ms) */
  double m2         = mc2mass2(vectorGetValue(parameter,"chirpmass"), eta);
  double mt         = m1 + m2;
  double dmm        = (m2-m1)/mt; /*  = (delta m) / mt  (dimensionless) */
  double log_mt     = log(mt) + log(Msun);  /* (in Kg) */
  double log_eta    = log(eta);
  double eta2       = eta * eta;
  double eta3       = eta2 * eta;
  double log_mu     = log_eta + log_mt;
  double log_omega0 = log(4.0*pi);
  double log_tau0   = log(1.0);
  double t, phi, psi;
  double taucoef = 80.37679297973246 + log_eta - log_mt;               /*  (4.17) or (11) */
  double log_tau, tau18, tau28, tau38, tau48, tau58, tau68, tau78;
  double ci  = cos(vectorGetValue(parameter,"inclination"));
  double ci2 = ci*ci,     ci4 = ci2*ci2,   ci6 = ci4*ci2;
  double si2 = (1.0-ci2), si  = sqrt(si2), si4 = si2*si2, si5 = si4*si;
  double h_plus, h_cross;
  double Hp00, Hp05, Hp10, Hp15, Hp20, Hp25;
  double Hc00, Hc05, Hc10, Hc15, Hc20, Hc25;
  double plus10a  = 0.1666666666666667*((19.0+9.0*ci2-2.0*ci4)-eta*(19.0-11.0*ci2-6.0*ci4));
  double plus10b  = 1.333333333333333*si2*(1.0+ci2)*(1.0-3.0*eta);
  double plus15a  = ((57.0 + 60.0*ci2-ci4) - 2.0*eta*(49.0-12.0*ci2-ci4));
  double plus15b  = 13.5*((73.0+40.0*ci2-9.0*ci4) - 2.0*eta*(25.0-8.0*ci2-9.0*ci4));
  double plus15c  = 312.5*(1.0-2.0*eta)*si2*(1.0+ci2);
  double plus20a  = 0.008333333333333333*((22.0+396.0*ci2+145.0*ci4-5.0*ci6)        
                    + 1.666666666666667*eta*(706.0-216.0*ci2-251.0*ci4+15.0*ci6)
	            -5.0*eta2*(98.0-108.0*ci2+7.0*ci4+5.0*ci6));
  double plus20b  = 0.1333333333333333*si2*((59.0+35.0*ci2-8.0*ci4)
	            -1.666666666666667*eta*(131.0+59.0*ci2-24.0*ci4)
	            +5.0*eta2*(21.0-3.0*ci2-8.0*ci4));
  double plus20c  = 2.025*(1.0-5.0*eta+5.0*eta2)*si4*(1.0+ci2);
  double plus20d  = (11.0+7.0*ci2+10.0*(5.0+ci2)*0.6931471805599453);
  double plus25a  = si*dmm*(0.3458984375-0.3255859375*ci2+0.02354600694444444*ci4-0.0001085069444444444*ci6
                            +eta*(2.66015625+0.01692708333333333*ci2-0.04557291666666666*ci4+0.0004340277777777778*ci6)
                            +eta2*(-0.3744574652777778+0.2190755208333333*ci2-0.0005425347222222222*ci4-0.0003255208333333333*ci6));
  double plus25b  = pi*(6.333333333333333+3.0*ci2-0.6666666666666666*ci4
                        +eta*(-5.333333333333333+4.666666666666667*ci2+2.0*ci4));
  double plus25c  = si*dmm*(3.4541015625-4.4876953125*ci2-2.9900390625*ci4+0.1423828125*ci6
                            +eta*(-18.61640625+4.31953125*ci2+6.05390625*ci4-0.56953125*ci6)
	                    +eta2*(5.6888671875-5.3255859375*ci2-0.3216796875*ci4+0.4271484375*ci6));
  double plus25d  = (-16.75516081914556*(1.0+ci2)*si2*(1.0-3.0*eta));
  double plus25e  = si*dmm*(-11.73231336805556+4.408094618055555*ci2+9.019639756944445*ci4-1.695421006944444*ci6
                            +eta*(31.73828125-17.63237847222222*ci2-20.88758680555556*ci4+6.781684027777778*ci6)
                            +eta2*(-12.95301649305556+13.22428385416667*ci2+4.814995659722222*ci4-5.086263020833333*ci6));
  double plus25f  = dmm*(2.553146701388889*si5*(1.0+ci2)*(1.0-4.0*eta+3.0*eta2));
  double plus25g  = (-1.8+2.8*ci2+1.4*ci4+eta*(19.2-1.6*ci2-5.6*ci4));
  double plus25h  = si2*(1.0+ci2)*(3.806430074027250-eta*17.58595688874842);
  double cross10a = (ci/3.0)*((17.0-4.0*ci2)-eta*(13.0-12.0*ci2));
  double cross10b = 2.666666666666667*(1.0-3.0*eta)*ci*si2;
  double cross15a = ((63.0-5.0*ci2)-2.0*eta*(23.0-5.0*ci2));
  double cross15b = 13.5*((67.0-15.0*ci2)-2.0*eta*(19.0-15.0*ci2));
  double cross15c = 312.5*(1.0-2.0*eta)*si2;
  double cross20a = (ci/60.0)*((68.0+226.0*ci2-15.0*ci4)+1.666666666666667*eta*(572.0-490.0*ci2+45.0*ci4)
                    -5.0*eta2*(56.0-70.0*ci2+15.0*ci4));
  double cross20b = 0.2666666666666667*ci*si2*((55.0-12.0*ci2)-1.666666666666667*eta*(119.0-36.0*ci2)
                    +5.0*eta2*(17.0-12.0*ci2));
  double cross20c = 4.05*(1.0-5.0*eta+5.0*eta2)*ci*si4;
  double cross25a = 1.2*si2*ci*eta;
  double cross25b = ci*(2.0-4.4*ci2+eta*(-30.8+18.8*ci2));
  double cross25c = ci*si2*(-7.612860148054500+eta*35.17191377749683);
  double cross25d = si*ci*dmm*(-0.1188802083333333+0.1641493055555556*ci2-0.001519097222222222*ci4
                               +eta*(3.033854166666667-0.4079861111111111*ci2+0.006076388888888889*ci4)
                               +eta2*(-0.2823350694444444+0.1306423611111111*ci2-0.004557291666666667*ci4));
  double cross25e = pi*ci*(11.33333333333333-2.666666666666667*ci2-eta*(6.666666666666667-8.0*ci2));
  double cross25f = si*ci*dmm*(4.883203125-9.42890625*ci2+0.664453125*ci4
                               +eta*(-30.5953125+24.440625*ci2-2.6578125*ci4)
                               +eta2*(7.383984375-8.90859375*ci2+1.993359375*ci4));
  double cross25g = si2*ci*(-33.51032163829112*(1.0-3.0*eta));
  double cross25h = dmm*si*ci*(-22.10828993055556+26.85546875*ci2-4.747178819444445*ci4
                               +eta*(58.05121527777778-77.03993055555556*ci2+18.98871527777778*ci4)
                               +eta2*(-21.83702256944444+36.07855902777778*ci2-14.24153645833333*ci4));
  double cross25i = dmm*si5*ci*(5.106293402777778*(1.0-4.0*eta+3.0*eta2));
  double sin1psi, sin2psi, sin3psi, sin4psi, sin5psi, sin6psi, sin7psi;
  double cos1psi, cos2psi, cos3psi, cos4psi, cos5psi, cos6psi, cos7psi;
  double constfactor = exp(-61.77448272506146 + log_mu - (vectorGetValue(parameter,"logdistance")+log(Mpc)));   /*  (6.01) */
  double x, sqrtx, oldx=0.0;
  double omega, omegacoef=exp(81.98623089216656 - log_mt); /*   =  (c^3)/(G*mt)  */
  double xcoef1 =  0.1842757936507937  + 0.2291666666666667*eta;
  double xcoef2 =  0.07709356890904510 + 0.1260799024470899*eta + 0.1076388888888889*eta2;
  double xcoef3 = -0.6948786875713585  + 0.1783508329381706*eta;
  double xcoef4 =  0.1189718784113057  + 2.584625427871787 *eta 
                   - 0.03438539858217592*eta2 + 0.07705500096450617*eta3;
  double xcoef5 = -0.825171564817328   - 0.6973257521615570*eta + 0.2393829775448565*eta2;
  double phicoef1 =  0.460689484126984  +  0.5729166666666666 *eta;
  double phicoef2 =  0.6418722070533943 +  1.103961278521825  *eta + 0.90576171875*eta2;
  double phicoef3 = -0.7057224708076262 +  0.09970875121256667*eta;
  double phicoef4 =  0.2268860596302182 - 19.53077295136832   *eta 
                     + 0.08423124040876116*eta2 - 0.6666536684389468*eta3;
  double phicoef5 =  3.415308237927678  +  2.975587931103962  *eta - 0.8629798504672995*eta2;
  double x_isco = 1.0/6.0; /* pow( (pi * f_isco)/omegacoef , 2.0/3.0); */
  int i, terminate=0;
  double *timedomainwaveform=NULL;
  double tc = vectorGetValue(parameter,"time");
  double phase = vectorGetValue(parameter,"phase");

  /* fill `timedomainwaveform' with time-domain template: */
  timedomainwaveform = (double*) malloc(sizeof(double)*DF->dataSize);
  for (i=0; i<DF->dataSize; ++i){
    /* determine time left until coalescence, "(t_c-t)" in (4.17)/(11): */
    t = (tc - DF->dataStart) - ((double)i)*DF->dataDeltaT; 
    if ((t>0.0) && (!terminate)) {  /*  (before t_c and before frequency reaches its maximum) */
      /*  determine `dimensionless time variable' tau: */
      log_tau = taucoef + log(t);                                                /*  (4.17), (11) */
      tau18   = exp(0.125 * log_tau);   /* = tau ^ (1/8) */
      tau28   = exp(0.25  * log_tau);   /* = tau ^ (2/8) */
      tau38   = exp(0.375 * log_tau);   /* = tau ^ (3/8) */
      tau48   = exp(0.5   * log_tau);   /* = tau ^ (4/8) */
      tau58   = exp(0.625 * log_tau);   /* = tau ^ (5/8) */
      tau68   = exp(0.75  * log_tau);   /* = tau ^ (6/8) */
      tau78   = exp(0.875 * log_tau);   /* = tau ^ (7/8) */
      /* determine (dimensionless) `frequency' x: */
      x = (0.25/tau28) * (1.0 + xcoef1/tau28 - 0.6283185307179586/tau38
                          + xcoef2/tau48 + xcoef3/tau58
                          + (xcoef4-0.03184523809523809*(log_tau-5.545177444479562))/tau68 
                          + xcoef5/tau78);                                        /*  (12)  */
      if ((x > x_isco) || (x < oldx)){  /* (frequency decreases  ==>  signal is terminated) */
        timedomainwaveform[i] = 0.0; 
        terminate = 1;
      }
      else {                    /*  (frequency still increasing  ==>  keep on computing...) */
        oldx    = x;
        sqrtx   = sqrt(x);
        /* derive angular frequency omega: (omega/pi gives frequency in Hz) */
        omega   = omegacoef*x*sqrtx;   /*  = ((c^3)/(G*mt)) * x^(3/2)                (4.13) */
        /* determine phase phi: */
	phi     = phase - (1.0/eta) * 
                  (tau58 + phicoef1*tau38 - 2.356194490192345*tau28
		   + phicoef2*tau18 + phicoef3*(log_tau-log_tau0)
                   + (phicoef4 + 0.2388392857142857*(log_tau-5.545177444479562))/tau18
                   + phicoef5/tau28);                                             /*  (13)    */
        /* derive `basic phase' psi: */
        /* psi     = phi - 2.0*x*sqrtx * (log(omega)-log_omega0); */              /*  (6.12)  */
	psi     = phi - 2.0*x*sqrtx * (log(omega)-log_omega0) * (1.0-(eta/2.0)*x); /* Arun et al. (5.6) */
	sin1psi = sin(psi);      cos1psi = cos(psi);
	sin2psi = sin(2.0*psi);  cos2psi = cos(2.0*psi);
	sin3psi = sin(3.0*psi);  cos3psi = cos(3.0*psi);
	sin4psi = sin(4.0*psi);  cos4psi = cos(4.0*psi);
	sin5psi = sin(5.0*psi);  cos5psi = cos(5.0*psi);
	sin6psi = sin(6.0*psi);  cos6psi = cos(6.0*psi);
	sin7psi = sin(7.0*psi);  cos7psi = cos(7.0*psi);
        /* determine PN plus- & cross-terms: */
	Hp00    = -(1.0+ci2)*cos2psi - (si2/96.0)*(17.0+ci2);                     /*  (6.02), Arun et al (5.7a) */
	Hp05    = -(si/8.0)*dmm * ((5.0+ci2)*cos1psi - 9.0*(1.0+ci2)*cos3psi);    /*  (6.03)  */
	Hp10    = plus10a*cos2psi - plus10b*cos4psi;                              /*  (6.04)  */
	Hp15    = (si/192.0)*dmm * (plus15a*cos1psi - plus15b*cos3psi + plus15c*cos5psi) 
                  - 6.283185307179586*(1.0+ci2)*cos2psi;                          /*  (6.05)  */
	Hp20    = plus20a*cos2psi + plus20b*cos4psi - plus20c*cos6psi
	          +si/40.0*dmm*(plus20d*sin1psi-15.70796326794897*(5.0+ci2)*cos1psi 
                  -79.52442081079560*(1.0+ci2)*sin3psi+424.115008234622*(1.0+ci2)*cos3psi);
                                                                                  /*  (6.06)  */
        Hp25    = cos1psi*plus25a + cos2psi*plus25b + cos3psi*plus25c
                  + cos4psi*plus25d + cos5psi*plus25e + cos7psi*plus25f
                  + sin2psi*plus25g + sin4psi*plus25h;                            /*  Arun & al. (5.09) */
	Hc00    = -2.0*ci*sin2psi;                                                /*  (6.07)  */
	Hc05    = -0.75*si*ci*dmm*(sin1psi-3.0*sin3psi);                          /*  (6.08)  */
	Hc10    = cross10a*sin2psi - cross10b*sin4psi;                            /*  (6.09)  */
	Hc15    = ((si*ci)/96.0)*dmm * 
                  (cross15a*sin1psi - cross15b*sin3psi + cross15c*sin5psi)
                  -12.56637061435917*ci*sin2psi;                                  /*  (6.10)  */
	Hc20    = cross20a*sin2psi + cross20b*sin4psi - cross20c*sin6psi
	          -0.15*si*ci*dmm*(9.931471805599454*cos1psi+15.70796326794897*sin1psi
	          -26.50814027026520*cos3psi - 141.3716694115407*sin3psi);        /*  (6.11)  */
        Hc25    = cross25a + cos2psi*cross25b + cos4psi*cross25c
                  + sin1psi*cross25d + sin2psi*cross25e + sin3psi*cross25f
                  + sin4psi*cross25g + sin5psi*cross25h + sin7psi*cross25i;       /*  Arun & al. (5.10) */
        /* and finally - the actual signal: */
	h_plus  = h_cross = constfactor * x;
	h_plus  *= Hp00 + sqrtx*(Hp05 + sqrtx*(Hp10 + sqrtx*(Hp15 + sqrtx*(Hp20 + sqrtx*Hp25))));
	h_cross *= Hc00 + sqrtx*(Hc05 + sqrtx*(Hc10 + sqrtx*(Hc15 + sqrtx*(Hc20 + sqrtx*Hc25))));/* (6.01) */
	timedomainwaveform[i] = Fplus*h_plus + Fcross*h_cross;                    /*  (3.10)  */
      }
    }
    else timedomainwaveform[i] = 0.0;  /*  (after t_c or after termination) */
  }
  /* FT'd waveform is returned in "output": */
  FTexec(DF, timedomainwaveform, output);
  free(timedomainwaveform);
}


void templateStatPhase(DataFramework *DF, vector *parameter, double Fplus, double Fcross,
                       double complex *output, double PNOrder)
/*************************************************************/
/* returns the (analytic) frequency-domain template.         */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* 2.5PN stationary-phase approximation                      */
/* following  Tanaka/Tagoshi (2000), Phys.Rev.D 62(8):082001 */
/* or Christensen/Meyer (2001), Phys.Rev.D 64(2):022001      */
/* 'PNOrder' may only be either ==2.0 or ==2.5.              */
/* Only difference between 2.0 and 2.5 PN order is that      */
/* for 2.0PN the 5th coefficient (a[4]) in Psi is left out.  */
/*************************************************************/
{
  double eta     = vectorGetValue(parameter,"massratio");
  double mt      = mc2mt(vectorGetValue(parameter,"chirpmass"), eta)*Msun;
  double log_q   = log(mt) + log(pi) + log(G) - 3.0*log(c);
  double log_eta = log(eta);
  double a[5];
  long i;
  double f, f01, f02, f04, f06, f07, f10, Psi, twopitc;
  double ampliConst;
  double sineCoef, cosineCoef;
  double complex cosinechirp;
  double phase = vectorGetValue(parameter,"phase");

  ampliConst = 0.5*log(5.0)+(5.0/6.0)*log(G)-log(2.0)-0.5*log(6.0)-(2.0/3.0)*log(pi)-1.5*log(c);
  ampliConst = exp(ampliConst+0.5*log_eta+(5.0/6.0)*log(mt)
                   -vectorGetValue(parameter,"logdistance")-log(Mpc));
  cosineCoef = Fplus  * (-0.5*(1.0+pow(cos(vectorGetValue(parameter,"inclination")),2.0)));
  sineCoef   = Fcross * (-1.0*cos(vectorGetValue(parameter,"inclination")));
  twopitc = 2.0 * pi * (vectorGetValue(parameter,"time") - DF->dataStart);
  a[0] =  exp(-3.75341797525151 - 1.66666666666667*log_q - log_eta);
  a[1] =  exp(log(44.2261904761905+55.0*eta) - 5.95064255258773 - log_eta - log_q);
  a[2] = -exp(0.163900632837675 - 0.666666666666667*log_q - log_eta);
  a[3] =  exp(-3.75341797525151 - log_eta - 0.333333333333333*log_q
              + log(30.1031529509952+53.859126984127*eta+42.8472222222222*exp(2.0*log_eta)));
  a[4] = (PNOrder<2.5) ? 0.0 : exp(-3.70730037807022-log_eta+log(153.353174603175+5.0*eta));

  /* fill the (complex-valued) "output[]" with template: */
  for (i=0; i<DF->FTSize; ++i){
    if ((i > DF->maxInd) || (i < DF->minInd)) /* (no computations outside freq. range) */
      output[i] = 0.0 + I*0.0;
    else {
      f    = ((double)i) * DF->FTDeltaF;
      f01  = pow(f, -1.0/6.0);             /* = f^-1/6  */
      f02  = f01*f01;                      /* = f^-2/6  */
      f04  = f02*f02;                      /* = f^-4/6  */
      f06  = f04*f02;                      /* = f^-6/6  */
      f07  = f06*f01;                      /* = f^-7/6  */
      f10  = f06*f04;                      /* = f^-10/6 */
      Psi = a[0]*f10 + a[1]*f06 + a[2]*f04 + a[3]*f02; /* + a[4]*log(f);  */
      if (PNOrder>2.0) /*  "a[4]*log(f)" - coefficient ignored for 2.0 PN order  */
        Psi += a[4]*log(f); 
      cosinechirp = ampliConst * f07 * cexp(-I*(Psi+phase+twopitc*f));
      output[i] = cosineCoef*cosinechirp + sineCoef*I*cosinechirp;
      output[i] /= DF->dataDeltaT;
    }
  }
}


void templateR2PN(DataFramework *DF, vector *parameter, double Fplus, double Fcross,
                  double complex *output)
/***************************************************************************/
/* returns the (numerically FT'd) frequency-domain template.               */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Restricted 2.0PN approximation as used for modelling MBH binaries       */
/* in the 1st & 2nd Mock LISA Data Challenge.                              */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Following the conventions of Section 3.4.6                              */
/* (pp.21 sqq) of `Document for Challenge 1',                              */
/* See: http://svn.sourceforge.net/viewvc/+checkout+/lisatools/Docs/challenge1.pdf */
/*      (replace "+" by "*" in above address)                              */
/* See also section 4.4 of                                                 */
/*    Arnaud et al. (2007): "An overview of the second                     */ 
/*    round of the Mock LISA Data Challenges".                             */
/*    Class. and Quantum Grav., 24(19):S551--S564.                         */
/* and references therein.                                                 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Note that this is (by now?) implemented without any tapering            */
/* (see equation (12) in Arnaud et al.).                                   */
/***************************************************************************/
{
  /*double MsunS = 4.92549095e-6;*/
  double t, logtau, Momega, Phi;
  double eta = vectorGetValue(parameter,"massratio");
  double m1 = mc2mass1(vectorGetValue(parameter,"chirpmass"), eta);
  double m2 = mc2mass2(vectorGetValue(parameter,"chirpmass"), eta);
  double totalmass = (m1+m2);              /* still in units if Msun */
  double reducedmass = eta*totalmass;      /* also in units if Msun  */
  double taufactor = 3.0*log(c)-log(G) + log(eta) - log(5.0) - log(totalmass) - log(Msun);
  double MomegaLSO  = 0.0138 * 1e6 * Msun;
  double MomegaMECO = (-((54.0+6.0*eta-6.0*sqrt(1539.0-(1008.0+19.0*eta)*eta))/(81.0-(57.0-eta)*eta))/27)*Msun;
  double maxMomega, oldMomega=-HUGE_VAL;
  double tau28, tau38, tau48;
  double omegacoef1 = (11.0/32.0)*eta + 743.0/2688.0;
  double omegacoef2 = 0.3*pi;
  double omegacoef3 = 1855099.0/14450688.0 + ((371.0/2048.0)*eta + 56975.0/258048.0)*eta;
  double Mo23, Mo43, Mo53;
  double Phicoef1 = 1.0/(32.0*eta);
  double Phicoef2 = 3715.0/1008.0 + (55.0/12.0)*eta;
  double Phicoef3 = 10.0*pi;
  double Phicoef4 = 15293365.0/1016064.0 + (27145.0/1008.0 + (3085.0/144.0)*eta)*eta;
  double hplus, hcross;
  double cosincli = cos(vectorGetValue(parameter,"inclination"));
  double hpcoef = exp(-61.77448272506146+log(reducedmass)+log(Msun)-vectorGetValue(parameter,"logdistance")-log(Mpc)) * (1.0+cosincli*cosincli);
  double hccoef = exp(-61.77448272506146+log(reducedmass)+log(Msun)-vectorGetValue(parameter,"logdistance")-log(Mpc)) * 2.0*cosincli;
  double tc = vectorGetValue(parameter,"time");
  double phase = vectorGetValue(parameter,"phase");
  int finished = 0;
  long i;
  double *timedomainwaveform=NULL;

  /* fill `timedomainwaveform' with time-domain template: */
  timedomainwaveform = (double*) malloc(sizeof(double)*DF->dataSize);
  maxMomega = MomegaLSO < MomegaMECO ? MomegaLSO : MomegaMECO;
  for (i=0; i<DF->dataSize; ++i){
    t = ((double) i) * DF->dataDeltaT;
    if ((t < tc) && (!finished)){
      logtau = taufactor + log((tc - DF->dataStart) - t);  /* (3.71) */
      tau28 = exp(0.25  * logtau);  /* tau^(2/8) */
      tau38 = exp(0.375 * logtau);  /* tau^(3/8) */
      tau48 = exp(0.5   * logtau);  /* tau^(4/8) */
      Momega = (0.125/tau38) * (1.0 + omegacoef1/tau28 - omegacoef2/tau38 + omegacoef3/tau48); /* (3.88) */
    }
    else Momega = HUGE_VAL;
    /* if before reaching LSO, derive waveform, otherwise set to zero: */
    if ((Momega < maxMomega) && (Momega>=oldMomega)){ 
      oldMomega = Momega;
      Mo23 = pow(Momega, 2.0/3.0);  /* Momega^(2/3) */
      Mo43 = Mo23 * Mo23;           /* Momega^(4/3) */
      Mo53 = Momega * Mo23;         /* Momega^(5/3) */
      Phi  = (-Phicoef1/Mo53) * (1.0 + Phicoef2*Mo23 - Phicoef3*Momega + Phicoef4*Mo43); /* (3.89) */
      Phi += phase;
      hplus  = hpcoef * Mo23 * cos(2.0*Phi);              /*  (3.92)  */
      hcross = hccoef * Mo23 * sin(2.0*Phi);              /*  (3.93)  */
    }
    else {
      finished = 1;
      hplus = hcross = 0.0;     /* (after reaching LSO, or after t_c) */
    }
    timedomainwaveform[i] = Fplus*hplus + Fcross*hcross;  /*  (3.10)  */
  }
  /* FT'd waveform is returned in "output": */
  FTexec(DF, timedomainwaveform, output);
  free(timedomainwaveform);
}



void templateLAL(DataFramework *DF, vector *parameter, double Fplus, double Fcross,
                 double complex *output, int approximant, int order)
/***************************************************************/
/* returns the (numerically FT'd) frequency-domain template.   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Using LAL inspiral templates.                               */
/* See also "LALInspiralTest.c".                               */
/***************************************************************/
{
  static LALStatus status;
  static InspiralTemplate params;
  static REAL4Vector *LALSignal=NULL;
  UINT4 n;
  double sineCoef, cosineCoef;
  long i;
  int FDomain;

  double mc  = vectorGetValue(parameter,"chirpmass"); 
  double eta = vectorGetValue(parameter,"massratio");
  double chirptime;
  double *timedomainwaveform=NULL;
  double timeshift = 0.0; /* time by which to shift template (in seconds) */
  double twopit;
  double complex *cosinechirp=NULL, *sinechirp=NULL;

  /* some (fixed) settings: */
  params.OmegaS      = 0.0;     /* (?) */
  params.Theta       = 0.0;     /* (?) */
  /* params.Zeta2    = 0.0; */  /* (?) */
  params.ieta        = 1; 
  params.nStartPad   = 0;
  params.nEndPad     = 0;
  params.massChoice  = m1Andm2;
  params.approximant = approximant;  /*  TaylorT1, ...              */
  params.order       = order;        /*  0=Newtonian, ..., 7=3.5PN  */
  params.fLower      = DF->minF * 0.9;
  params.fCutoff     = (DF->FTSize-1)*DF->FTDeltaF;  /* (Nyquist freq.) */
  params.tSampling   = 1.0/DF->dataDeltaT;
  params.startTime   = 0.0;

  /* actual inspiral parameters: */
  params.mass1       = mc2mass1(mc, eta);
  params.mass2       = mc2mass2(mc, eta);
  params.startPhase  = vectorGetValue(parameter,"phase");
  /* set distance according to approximant's specific requirements (ARRR): */
  if ((params.approximant == EOB) 
      || (params.approximant == EOBNR)
    /*|| (params.approximant == BBHPhenTD)
      || (params.approximant == BBHPhenFD)*/
      || (params.approximant == TaylorT3)
      || (params.approximant == IMRPhenomA)
      )
    params.distance  = exp(vectorGetValue(parameter,"logdistance")+log(Mpc)); /* distance in metres */
  else if ((params.approximant == TaylorT1)
           || (params.approximant == TaylorT2)
           || (params.approximant == PadeT1)
           || (params.approximant == TaylorF1)
           || (params.approximant == TaylorF2)
           || (params.approximant == PadeF1)
           || (params.approximant == BCV)
	   )
    params.distance  = exp(vectorGetValue(parameter,"logdistance"));          /* distance in Mpc */
  else                                                     
    params.distance  = exp(vectorGetValue(parameter,"logdistance")+log(Mpc)-log(c)); /* distance in seconds */

  cosineCoef = Fplus  * (-0.5*(1.0+pow(cos(vectorGetValue(parameter,"inclination")),2.0)));
  sineCoef   = Fcross * (-1.0*cos(vectorGetValue(parameter,"inclination")));

  /* ensure proper "fCutoff" setting: */
  if (params.fCutoff >= 0.5*params.tSampling)
    params.fCutoff = 0.5*params.tSampling - 0.5*DF->FTDeltaF;
  if (! (params.tSampling > 2.0*params.fCutoff)){
    printf(" : WARNING: 'LALInspiralSetup()' (called within 'LALInspiralWavelength()')\n");
    printf(" :          requires (tSampling > 2 x fCutoff) !!\n");
  }

  /* ensure compatible sampling rate: */
  if ((params.approximant == EOBNR)
      && (fmod(log((double)params.tSampling)/log(2.0),1.0) != 0.0)) {
    printf(" : WARNING: \"EOBNR\" templates require power-of-two sampling rates!!\n");
    printf(" :          (params.tSampling = %f Hz)\n", params.tSampling);
  }

  /* compute other elements of `params', check out the `.tC' value, */
  /* shift the start time to match the coalescence time,            */
  /* and eventually re-do parameter calculations:                   */

  /*printf(":: LALInspiralWave(..., approximant=%d, order=%d)\n", params.approximant, params.order);
    printf(" :  tC way before: %f\n", params.tC);*/


  LALInspiralParameterCalc(&status, &params);
  /*printf(" :  tC wee before: %f\n", params.tC);*/
  chirptime = params.tC;
  if ((params.approximant != TaylorF2) && (params.approximant != BCV)) {
    params.startTime = (vectorGetValue(parameter,"time") - DF->dataStart) - chirptime;
    LALInspiralParameterCalc(&status, &params); /* (re-calculation necessary? probably not...) */
  }

  /* compute "params.signalAmplitude" slot: */
  LALInspiralRestrictedAmplitude(&status, &params);

  /* figure out inspiral length & set `n': */
  /* LALInspiralWaveLength(&status, &n, params); */
  n = DF->dataSize;
  /* allocate waveform vector: */
  LALCreateVector(&status, &LALSignal, n);
for (i=0; i<DF->dataSize; ++i) LALSignal->data[i] = 0.0;
  /* compute actual waveform: */

/*rams.tC = 0.0;
 printf(" :  tC before    : %f\n", params.tC);*/

  /* REPORTSTATUS(&status); */
  LALInspiralWave(&status, LALSignal, &params);
  /* REPORTSTATUS(&status); */

  /*printf(" :  tC after     : %f\n", params.tC);/*
  /* REPORTSTATUS(&status);*/
  /* frequency domain or time domain waveform? */
  FDomain = ((params.approximant == TaylorF1)
             || (params.approximant == TaylorF2)
             || (params.approximant == PadeF1)
             || (params.approximant == BCV)
             /*|| (params.approximant == BBHPhenFD)*/
             );

  if (FDomain && (n % 2 != 0))
    printf(" : WARNING: frequency-domain LAL waveforms require even number of samples `n' !!\n");
             
  cosinechirp = (double complex*) malloc(sizeof(double complex)*DF->FTSize);
  if (!FDomain){ /* (approximant yields a time-domain waveform) */
    /* fill `timedomainwaveform' with time-domain (cosine chirp) waveform: */
    timedomainwaveform = (double*) malloc(sizeof(double)*DF->dataSize);
    for (i=0; i<DF->dataSize; ++i)
      /* timedomainwaveform[i] = (i<n) ? (LALSignal->data[i]) : 0.0; */
      timedomainwaveform[i] = LALSignal->data[i];
    LALDestroyVector(&status, &LALSignal);

    /* numerically Fourier-transform waveform: */
    FTexec(DF, timedomainwaveform, cosinechirp);
    free(timedomainwaveform);
  }
  else { /* (approximant yields a frequency-domain waveform) */
    /* copy over: */
    cosinechirp[0] = LALSignal->data[0];
    for (i=1; i<DF->FTSize-1; ++i)
      cosinechirp[i] = (LALSignal->data[i] + I*LALSignal->data[DF->dataSize-i]);
    cosinechirp[DF->FTSize-1] = LALSignal->data[DF->FTSize-1];
    LALDestroyVector(&status, &LALSignal);
    /* normalise: */
    for (i=0; i<DF->FTSize; ++i)
      cosinechirp[i] *= DF->dataSize;  /*   * DF->dataDeltaT;   */
  }


  /* figure out TIME SHIFT (if necessary): */
  if ((params.approximant == TaylorT2) 
      || (params.approximant == TaylorF2))
    timeshift = (vectorGetValue(parameter,"time") - DF->dataStart) - chirptime;
  else if (params.approximant == BCV)
    timeshift = (vectorGetValue(parameter,"time") - DF->dataStart) - (((double)DF->dataSize)*DF->dataDeltaT);
  else if (params.approximant == IMRPhenomA)
    timeshift = 0.0; /*(vectorGetValue(parameter,"time") - DF->dataStart) - params.tC;*/

  /* time-shift the template: */
  if (timeshift != 0.0) { 
    twopit = 2.0*pi*timeshift;
    for (i=1; i<DF->FTSize-1; ++i)
      cosinechirp[i] *= cexp(-twopit*I*(((double)i) * DF->FTDeltaF));
  }

  /* determine eventual actual response (returned in "output"): */
  for (i=0; i<DF->FTSize; ++i)
    output[i] = cosineCoef*cosinechirp[i] + sineCoef*I*cosinechirp[i];
  free(cosinechirp);

  /*
   * NOTE: the dirty trick here is to assume the LAL waveform to constitute
   *       the cosine chirp and then derive the corresponding sine chirp 
   *       as the orthogonal ("i x cosinechirp") waveform.
   *       In general they should not necessarily be only related 
   *       by a mere phase shift though...
   */
}



void templateSineGaussianBurst(DataFramework *DF, vector *parameter, double Fplus, double Fcross,
                               double complex *output)
/*************************************************************/
/* returns the (numerically FT'd) frequency-domain template. */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Sine-Gaussian burst signal.                               */
/* Parameters (8):                                           */
/*   - azimuth                                               */
/*   - altitude                                              */
/*   - polarisation                                          */
/*   - time                                                  */
/*   - logamplitude                                          */
/*   - logsigma                                              */
/*   - frequency                                             */
/*   - phase                                                 */
/*************************************************************/
{
  double hplus, hcross;
  double t, tsigma;
  double f = vectorGetValue(parameter,"frequency");
  double mu = vectorGetValue(parameter,"time");
  double sigma = exp(vectorGetValue(parameter,"logsigma"));
  double logampli = vectorGetValue(parameter,"logamplitude");
  double phase = vectorGetValue(parameter,"phase");
  double twopif = 2.0*pi*f;
  long i;
  double *timedomainwaveform=NULL;

  /* fill `timedomainwaveform' with time-domain template: */
  timedomainwaveform = (double*) malloc(sizeof(double)*DF->dataSize);

  for (i=0; i<DF->dataSize; ++i){
    t = ((double)i)*DF->dataDeltaT - (mu - DF->dataStart);  /* (time from t_c) = "(t-t_c)" */
    /* ==  (tstart + i*deltaT) - t_c */
    tsigma = t / sigma; /* = t in units of sigma */
    if (fabs(tsigma)<39.0) { /* otherwise set template to zero */
      hplus = hcross = exp(logampli - 0.5*tsigma*tsigma);
      if (f > 0.0) {
        hplus  *= cos(twopif*t - phase);
        hcross *= 0.0; /* (one polarisation only) */
      }
      timedomainwaveform[i] = Fplus*hplus + Fcross*hcross;
    }
    else timedomainwaveform[i] = 0.0; 
  }

  /* FT'd waveform is returned in "output": */
  FTexec(DF, timedomainwaveform, output);
  free(timedomainwaveform);
}


void inject(DataFramework *DF, int coherentN, int waveform, vector *parameter)
     /* Inject a signal into the (already FT'd) data. */
{
  double complex *FourierTemplate=NULL;
  double *indivSNR;
  double netSNR;
  long i, j;
  double locdeltat, locpolar, locazi, localti;
  vector localparameter;

  vectorInit(&localparameter);
  for (j=0; j<parameter->dimension; ++j)
    if ((strcmp(parameter->name[j],"longitude")!=0) 
        && (strcmp(parameter->name[j],"latitude")!=0))
      vectorAdd(&localparameter, parameter->name[j], parameter->value[j]);
  vectorAdd(&localparameter, "azimuth", 0.0);
  vectorAdd(&localparameter, "altitude", 0.0);

  for (i=0; i<coherentN; ++i) {
    FourierTemplate = (double complex*) malloc(sizeof(double complex)*DF[i].FTSize);
    /* determine "local" parameters: */
    localParameters(parameter, DF[i].ifo, &locdeltat, &locpolar, &localti, &locazi);
    vectorSetValue(&localparameter, "time",         vectorGetValue(parameter,"time")+locdeltat);
    vectorSetValue(&localparameter, "polarisation", locpolar);
    vectorSetValue(&localparameter, "azimuth",      locazi);
    vectorSetValue(&localparameter, "altitude",     localti);

    /* compute Fourier-domain template: */
    signaltemplate(&DF[i], waveform, &localparameter, FourierTemplate);
    /* compute sum-of-squares: */
    for (j=0; j<=DF[i].FTSize; ++j)
      DF[i].dataFT[j] += FourierTemplate[j];
    free(FourierTemplate);
  }
  if (verbose) {
    printf(" | injected '#%d' waveform signal.\n", waveform);
    indivSNR = (double*) malloc(sizeof(double)*coherentN);
    netSNR = signaltonoiseratio(DF, coherentN, waveform, parameter, indivSNR);
    printf(" | SNR: %.3f (%.3f", netSNR, indivSNR[0]);
    if (coherentN>1) for (i=1; i<coherentN; ++i) printf(", %.3f", indivSNR[i]);
    printf(")\n");
    free(indivSNR);
  }
}


void dumptemplates(DataFramework *DF, vector *parameter, char *filenameF, char *filenameT)
/* Write the different signal templates                 */
/* and (1-sided) power spectral density to a text file. */
/* (for sanity checking etc.)                           */
{
  double complex *FourierTemplate01=NULL;
  double complex *FourierTemplate02=NULL;
  double complex *FourierTemplate03=NULL;
  double complex *FourierTemplate04=NULL;
  double complex *FourierTemplate05=NULL;
  double complex *FourierTemplate06=NULL;

  fftw_plan InvFTplan;
  double complex *fourierdomain=NULL;
  double *timedomain=NULL;

  double *TimeTemplate01=NULL;
  double *TimeTemplate02=NULL;
  double *TimeTemplate03=NULL;
  double *TimeTemplate04=NULL;
  double *TimeTemplate05=NULL;
  double *TimeTemplate06=NULL;
  long i;
  FILE *textfile;

  double locdeltat, locpolar, locazi, localti;
  vector localparameter;
  vectorInit(&localparameter);
  for (i=0; i<parameter->dimension; ++i)
    if ((strcmp(parameter->name[i],"longitude")!=0) 
        && (strcmp(parameter->name[i],"latitude")!=0))
      vectorAdd(&localparameter, parameter->name[i], parameter->value[i]);
  vectorAdd(&localparameter, "azimuth", 0.0);
  vectorAdd(&localparameter, "altitude", 0.0);
  localParameters(parameter, DF[0].ifo, &locdeltat, &locpolar, &localti, &locazi);
  vectorSetValue(&localparameter, "time",         vectorGetValue(parameter,"time")+locdeltat);
  vectorSetValue(&localparameter, "polarisation", locpolar);
  vectorSetValue(&localparameter, "azimuth",      locazi);
  vectorSetValue(&localparameter, "altitude",     localti);

  printf(" : writing F'domain templates to file '%s'...\n", filenameF);
  FourierTemplate01 = (double complex*) malloc(sizeof(double complex)*DF->FTSize);
  FourierTemplate02 = (double complex*) malloc(sizeof(double complex)*DF->FTSize);
  FourierTemplate03 = (double complex*) malloc(sizeof(double complex)*DF->FTSize);
  FourierTemplate04 = (double complex*) malloc(sizeof(double complex)*DF->FTSize);
  FourierTemplate05 = (double complex*) malloc(sizeof(double complex)*DF->FTSize);
  FourierTemplate06 = (double complex*) malloc(sizeof(double complex)*DF->FTSize);

  /* compute (Fourier-domain) templates: */
  signaltemplate(DF, i25SP,          &localparameter, FourierTemplate01);
  signaltemplate(DF, i2535,          &localparameter, FourierTemplate02);
  signaltemplate(DF, iR2PN,          &localparameter, FourierTemplate03);
  signaltemplate(DF, iLALTT2PN20,    &localparameter, FourierTemplate04);
  signaltemplate(DF, iLALTT3PN20,    &localparameter, FourierTemplate05);
  signaltemplate(DF, iLALIMRPhenomA, &localparameter, FourierTemplate06);

  /* write to file: */
  textfile = fopen(filenameF, "w");
  fprintf(textfile,"f PSD real01 imag01 real02 imag02 real03 imag03 real04 imag04 real05 imag05 real06 imag06\n");
  for (i=0; i<DF->FTSize; ++i){
    fprintf(textfile, "%f %e", ((double)i)*DF->FTDeltaF, exp(DF->powspec[i]));
    fprintf(textfile, " %e %e", creal(FourierTemplate01[i]), cimag(FourierTemplate01[i]));
    fprintf(textfile, " %e %e", creal(FourierTemplate02[i]), cimag(FourierTemplate02[i]));
    fprintf(textfile, " %e %e", creal(FourierTemplate03[i]), cimag(FourierTemplate03[i]));
    fprintf(textfile, " %e %e", creal(FourierTemplate04[i]), cimag(FourierTemplate04[i]));
    fprintf(textfile, " %e %e", creal(FourierTemplate05[i]), cimag(FourierTemplate05[i]));
    fprintf(textfile, " %e %e", creal(FourierTemplate06[i]), cimag(FourierTemplate06[i]));
    fprintf(textfile,"\n");
  }
  fclose(textfile);
  printf(" : ...done.\n");

  printf(" : back-transforming Fourier-domain templates to time domain...\n");
  fourierdomain = (double complex*) malloc(sizeof(double complex)*DF->FTSize);
  timedomain = (double*) malloc(sizeof(double)*DF->dataSize);
  InvFTplan = fftw_plan_dft_c2r_1d(DF->dataSize, fourierdomain, timedomain, FFTW_ESTIMATE);

  for (i=0; i<DF->FTSize; ++i) fourierdomain[i] = FourierTemplate01[i];
  fftw_execute(InvFTplan);
  free(FourierTemplate01);
  TimeTemplate01 = (double*) malloc(sizeof(double)*DF->dataSize);
  for (i=0; i<DF->dataSize; ++i) TimeTemplate01[i] = timedomain[i]/DF->dataSize;

  for (i=0; i<DF->FTSize; ++i) fourierdomain[i] = FourierTemplate02[i];
  fftw_execute(InvFTplan);
  free(FourierTemplate02);
  TimeTemplate02 = (double*) malloc(sizeof(double)*DF->dataSize);
  for (i=0; i<DF->dataSize; ++i) TimeTemplate02[i] = timedomain[i]/DF->dataSize;

  for (i=0; i<DF->FTSize; ++i) fourierdomain[i] = FourierTemplate03[i];
  fftw_execute(InvFTplan);
  free(FourierTemplate03);
  TimeTemplate03 = (double*) malloc(sizeof(double)*DF->dataSize);
  for (i=0; i<DF->dataSize; ++i) TimeTemplate03[i] = timedomain[i]/DF->dataSize;

  for (i=0; i<DF->FTSize; ++i) fourierdomain[i] = FourierTemplate04[i];
  fftw_execute(InvFTplan);
  free(FourierTemplate04);
  TimeTemplate04 = (double*) malloc(sizeof(double)*DF->dataSize);
  for (i=0; i<DF->dataSize; ++i) TimeTemplate04[i] = timedomain[i]/DF->dataSize;

  for (i=0; i<DF->FTSize; ++i) fourierdomain[i] = FourierTemplate05[i];
  fftw_execute(InvFTplan);
  free(FourierTemplate05);
  TimeTemplate05 = (double*) malloc(sizeof(double)*DF->dataSize);
  for (i=0; i<DF->dataSize; ++i) TimeTemplate05[i] = timedomain[i]/DF->dataSize;

  for (i=0; i<DF->FTSize; ++i) fourierdomain[i] = FourierTemplate06[i];
  fftw_execute(InvFTplan);
  free(FourierTemplate06);
  TimeTemplate06 = (double*) malloc(sizeof(double)*DF->dataSize);
  for (i=0; i<DF->dataSize; ++i) TimeTemplate06[i] = timedomain[i]/DF->dataSize;

  free(InvFTplan);  free(fourierdomain);  free(timedomain);
  printf(" : ...done.\n");

  printf(" : writing T'domain templates to '%s'...\n", filenameT);
  textfile = fopen(filenameT, "w");
  fprintf(textfile,"t h01 h02 h03 h04 h05 h06\n");
  for (i=0; i<DF->dataSize; ++i){
    fprintf(textfile, "%f",  DF->dataStart + ((double)i)*DF->dataDeltaT);
    fprintf(textfile, " %e", TimeTemplate01[i]);
    fprintf(textfile, " %e", TimeTemplate02[i]);
    fprintf(textfile, " %e", TimeTemplate03[i]);
    fprintf(textfile, " %e", TimeTemplate04[i]);
    fprintf(textfile, " %e", TimeTemplate05[i]);
    fprintf(textfile, " %e", TimeTemplate06[i]);
    fprintf(textfile,"\n");
  }
  fclose(textfile);

  free(TimeTemplate01);
  free(TimeTemplate02);
  free(TimeTemplate03);
  free(TimeTemplate04);
  free(TimeTemplate05);
  free(TimeTemplate06);
  printf(" : ...done.\n");
}


double loglikelihoodOLD(DataFramework *DF, int waveform, vector *parameter)
/* Generic function to compute the likelihood.                 */
/* Uses the general "signaltemplate()" function to compute the */
/* frequency-domain waveform template corresponding to the     */
/* waveform model specified by the "waveform" argument, and    */
/* then matches data & signal using the power spectrum         */
/* and frequency range settings defined in the "DF" structure. */
{
  double complex *FourierTemplate=NULL;
  long i;
  double absdiff, chisquared=0.0;
  double logfactor = log(DF->dataDeltaT) - log(DF->dataSize) + log(2.0);
  FourierTemplate = (double complex*) malloc(sizeof(double complex)*DF->FTSize);
  /* compute Fourier-domain template: */
  signaltemplate(DF, waveform, parameter, FourierTemplate);
  /* compute sum-of-squares: */
  for (i=DF->minInd; i<=DF->maxInd; ++i){
    absdiff    =  cabs(DF->dataFT[i] - FourierTemplate[i]);
    chisquared += exp(logfactor + 2.0*log(absdiff) - DF->powspec[i]);
  }
  free(FourierTemplate);
  /* log-likelihood is negative sum-of-squares: */
  return -1.0 * chisquared;
}


double loglikelihood(DataFramework *DF, int coherentN, int waveform, vector *parameter)
/* Generic function to compute the likelihood.                 */
/* Uses the general "signaltemplate()" function to compute the */
/* frequency-domain waveform template corresponding to the     */
/* waveform model specified by the "waveform" argument, and    */
/* then matches data & signal using the power spectrum         */
/* and frequency range settings defined in the "DF" structure. */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Supplying "NULL" for the "parameter" will make the function */
/* NOT compute a template but rather plug in zeroes and com-   */
/* pute the "null"-likelihood instead.                         */
{
  double complex *FourierTemplate=NULL;
  long i, j, maxftsize;
  double absdiff, chisquared=0.0;
  double logfactor;
  double locdeltat, locpolar, locazi, localti;
  vector localparameter;
  maxftsize = 0;
  for (i=0; i<coherentN; ++i)
    if (DF[i].FTSize > maxftsize) maxftsize = DF[i].FTSize;
  FourierTemplate = (double complex*) malloc(sizeof(double complex)*maxftsize);

  vectorInit(&localparameter);
  if (parameter != NULL) {
    /* copy over everything except longitude & latitude, and add azimuth & altitude instead: */
    for (j=0; j<parameter->dimension; ++j)
      if ((strcmp(parameter->name[j],"longitude")!=0) 
          && (strcmp(parameter->name[j],"latitude")!=0))
        vectorAdd(&localparameter, parameter->name[j], parameter->value[j]);
    vectorAdd(&localparameter, "azimuth", 0.0);
    vectorAdd(&localparameter, "altitude", 0.0);
  }
  else { /* compute "NULL"-likelihood, assuming no signal present: */
    for (j=0; j<maxftsize; ++j) {
      FourierTemplate[j] = 0.0 + 0.0*I;
    }
  }
  /* loop over individual data sets / interferometers: */
  for (i=0; i<coherentN; ++i){
    logfactor = log(DF[i].dataDeltaT) - log(DF[i].dataSize) + log(2.0);
    if (parameter != NULL) {
      /* determine "local" parameters:    */
      localParameters(parameter, DF[i].ifo, &locdeltat, &locpolar, &localti, &locazi);
      vectorSetValue(&localparameter, "time",         vectorGetValue(parameter,"time")+locdeltat);
      vectorSetValue(&localparameter, "polarisation", locpolar);
      vectorSetValue(&localparameter, "azimuth",      locazi);
      vectorSetValue(&localparameter, "altitude",     localti);
      /* compute Fourier-domain template: */
      signaltemplate(&DF[i], waveform, &localparameter, FourierTemplate);
    }

    /* compute sum-of-squares:          */
    for (j=DF[i].minInd; j<=DF[i].maxInd; ++j){
      absdiff    =  cabs(DF[i].dataFT[j] - FourierTemplate[j]);
      chisquared += exp(logfactor + 2.0*log(absdiff) - DF[i].powspec[j]);
    }
  }

  free(FourierTemplate);
  vectorDispose(&localparameter);

  /* log-likelihood is negative sum-of-squares: */
  return -1.0 * chisquared;
}


double signaltonoiseratio(DataFramework *DF, int coherentN, int waveform, vector *parameter,
                          double indivSNRs[])
/* Computes signal-to-noise ratio (SNR).                                  */
/* If "indivSNRs" is not NULL, individual SNRs are returned here as well. */
{
  double complex *FourierTemplate=NULL;
  long i, j, maxftsize;
  double power, sum=0.0, networksum=0.0;
  double snrcoef;
  double locdeltat, locpolar, locazi, localti;
  vector localparameter;

  maxftsize = 0;
  for (i=0; i<coherentN; ++i)
    if (DF[i].FTSize > maxftsize) maxftsize = DF[i].FTSize;
  FourierTemplate = (double complex*) malloc(sizeof(double complex)*maxftsize);

  /* copy over everything except longitude & latitude, and add azimuth & altitude instead: */
  vectorInit(&localparameter);
  for (j=0; j<parameter->dimension; ++j)
    if ((strcmp(parameter->name[j],"longitude")!=0) 
        && (strcmp(parameter->name[j],"latitude")!=0))
      vectorAdd(&localparameter, parameter->name[j], parameter->value[j]);
  vectorAdd(&localparameter, "azimuth", 0.0);
  vectorAdd(&localparameter, "altitude", 0.0);

  /* loop over individual data sets / interferometers: */
  networksum = 0.0;
  for (i=0; i<coherentN; ++i){
    sum = 0.0;
    /* determine "local" parameters: */
    localParameters(parameter, DF[i].ifo, &locdeltat, &locpolar, &localti, &locazi);
    vectorSetValue(&localparameter, "time",         vectorGetValue(parameter,"time")+locdeltat);
    vectorSetValue(&localparameter, "polarisation", locpolar);
    vectorSetValue(&localparameter, "azimuth",      locazi);
    vectorSetValue(&localparameter, "altitude",     localti);

    /* compute Fourier-domain template: */
    signaltemplate(&DF[i], waveform, &localparameter, FourierTemplate);

    snrcoef = log(DF[i].dataDeltaT) - log((double)DF[i].dataSize);
    for (j=DF[i].minInd; j<=DF[i].maxInd; ++j){
      power      = 2.0 * log(cabs(FourierTemplate[j]));
      power      += snrcoef;
      sum        += exp(power - DF[i].powspec[j]);
    }
    if (indivSNRs != NULL)
      indivSNRs[i] = sqrt(4.0 * sum);
    networksum += sum;
  }

  free(FourierTemplate);
  vectorDispose(&localparameter);

  return sqrt(4.0 * networksum);
}


void clearMF(McmcFramework *MF)
/* cleans up pointers etc. in a "McmcFramework" structure. */
{
  if (verbose) printf(" | cleaning up 'McmcFramework' structure ...");
  vectorDispose(&MF->fixed);
  vectorDispose(&MF->startvalue);
  vectorDispose(&MF->priorparameters);
  free(MF->covariance);
  gsl_matrix_free(MF->covcholesky);
  free(MF);
  MF = NULL;
  if (verbose) printf(" done.\n");
}


int parindex(McmcFramework *MF, char parametername[])
/* returns index (w.r.t. that MF->startvalue vector) */
/* of parameter named 'parametername'.               */
{
  int i, index=-1;
  int notfound=1;
  if (MF->startvalue.dimension > 0) {
    i=0;
    while ((i<MF->startvalue.dimension) 
           && (notfound=(strcmp(MF->startvalue.name[i], parametername)!=0))) 
      ++i;
    if (notfound)
      printf(" : ERROR: requesting unknown vector element in 'parindex(...,\"%s\")'!\n", 
             parametername);
    else
      index = i;
  }
  else
    printf(" : ERROR: empty 'startvalue' vector in 'parindex(...,\"%s\")'!\n", 
           parametername);
  return index;
}

void setcov(McmcFramework *MF, int parameter1, int parameter2, double covariance)
/* Set an element of the proposal covariance matrix. */
{
  int pardim = MF->startvalue.dimension;
  if ((parameter1>=pardim) | (parameter2>=pardim))
    printf(" : ERROR: attempt to access invalid matrix element (%d,%d) in 'setcov()'!\n", 
           parameter1, parameter2);
  MF->covariance[parameter1 + parameter2*pardim] = covariance;
  MF->covariance[parameter2 + parameter1*pardim] = covariance;
}

void setvar(McmcFramework *MF, int parameter, double variance)
/* Set a main diagonal element of the proposal covariance matrix. */
{
  if (variance<=0.0)
    printf(" : ERROR: invalid variance value (%e) in 'setvar()'!\n", variance);
  setcov(MF, parameter, parameter, variance);
}

void setstdev(McmcFramework *MF, int parameter, double standarddeviation)
/* Set a main diagonal element of the proposal covariance matrix */
/* in terms of standard deviation.                               */
{
  setvar(MF, parameter, standarddeviation*standarddeviation);  
}

double getcov(McmcFramework *MF, int parameter1, int parameter2)
/* Retrieve an element of the proposal covariance matrix. */
{
  return MF->covariance[parameter1 + parameter2*MF->startvalue.dimension];
}

double getvar(McmcFramework *MF, int parameter)
/* Retrieve a main diagonal element of the proposal covariance matrix. */
{
  return getcov(MF, parameter, parameter);
}

double getstdev(McmcFramework *MF, int parameter)
/* Retrieve a main diagonal element of the proposal covariance matrix */
/* in terms of standard deviation.                                    */
{
  return sqrt(getvar(MF, parameter));
}

void setcor(McmcFramework *MF, int parameter1, int parameter2, double correlation)
/* Set an off-diagonal element of the proposal covariance matrix         */
/* in terms of the correlation coefficient.                              */
/* Requires the corresponding main diagonal elements to be set already!! */
{
  double covariance;
  if (fabs(correlation) > 1.0)
    printf(" : ERROR: Invalid correlation value (%e) in 'setcor()'!\n",correlation);
  if ((getstdev(MF, parameter1)<=0.0) | (getstdev(MF, parameter2)<=0.0))
    printf(" : WARNING: zero variance implies zero covariance in 'setcor(...,%d,%d,...)'!\n",
           parameter1, parameter2);
  covariance = getstdev(MF, parameter1) * getstdev(MF, parameter2) * correlation;
  setcov(MF, parameter1, parameter2, covariance);
}


void RandMVNorm(McmcFramework *MF, double *result)
/* generate multivariate-normal distributed random draw           */
/* with zero mean and covariance determined by `MF->covariance'   */
/* (or actually its Cholesky decomposition `MF->covcholesky',     */
/* which was initialised in the `init()' function).               */
/* The result is returned in the (pre-allocated!) `result' vector */
/* of dimension `MF->pardim'.                                     */
{
  int i;
  gsl_vector *draw;
  draw = gsl_vector_alloc(MF->pardim);
  /* first generate vector of i.i.d. standard normal draws: */
  for (i=0; i<MF->pardim; ++i)
    gsl_vector_set(draw, i, gsl_ran_ugaussian(GSLrandom));
  /* matrix-vector product: */
  gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, MF->covcholesky, draw);
  /* copy result to output: */
  for (i=0; i<MF->pardim; ++i)
    result[i] = gsl_vector_get(draw, i);
}


double logprior(McmcFramework *MF, vector *parameter)
/* 'wrapper' function for prior density.                    */
/* Checks the "MF->parameterset" character string and calls */
/* the corresponding log-prior-density function.            */
{
  double result = -HUGE_VAL;
  if (MF->parameterset == InspiralNoSpin)
    result = logpriorInspiralNospin(MF, parameter);
  else if (MF->parameterset == BurstSineGaussian)
    result = logpriorSineGaussianBurst(MF, parameter);
  else
    printf(" : ERROR: referring to undefined parameterset in 'logprior()'!\n");
  return(result);
}


double logpriorInspiralNospin(McmcFramework *MF, vector *parameter)
/* log - prior density for inspiral signal w/ no spin (9 parameters).              */
/* (not normalised!)                                                               */
/* Accounts for (approximated) detection probability based on intrinsic parameters */
/* as described in Roever (2007), http://hdl.handle.net/2292/2356, page 78 & sqq.  */
{
  double ampli90 = log(0.5) + (5.0/6.0)*log(4.0) + 0.5*log(8.0) 
                   - log(vectorGetValue(&MF->priorparameters,"dist90"));
  double ampli10 = log(0.5) + (5.0/6.0)*log(4.0) + 0.5*log(8.0) 
                   - log(vectorGetValue(&MF->priorparameters,"dist10"));
  double a = (ampli10+ampli90)/2.0;
  double b = (ampli10-ampli90)/(-2.0*log(0.1/0.9));
  double ampli;  /* logarithmic amplitude */
  double prior=0.0;
  double cosiota = cos(vectorGetValue(parameter,"inclination"));
  double cosiota2 = cosiota*cosiota;

  /* check whether parameter values fall within domain: */
  if ((vectorGetValue(parameter,"massratio")>0.25) | (vectorGetValue(parameter,"massratio")<=0.0)
      | (mc2mass1(vectorGetValue(parameter,"chirpmass"), vectorGetValue(parameter,"massratio")) < vectorGetValue(&MF->priorparameters,"massLower"))
      | (mc2mass2(vectorGetValue(parameter,"chirpmass"), vectorGetValue(parameter,"massratio")) > vectorGetValue(&MF->priorparameters,"massUpper"))
      | (vectorGetValue(parameter,"time") < vectorGetValue(&MF->priorparameters,"timeLower")) 
      | (vectorGetValue(parameter,"time") > vectorGetValue(&MF->priorparameters,"timeUpper"))
      | (vectorGetValue(parameter,"longitude") < -pi) | (vectorGetValue(parameter,"longitude") > pi)
      | (fabs(vectorGetValue(parameter,"latitude")) > pi/2.0)
      | (vectorGetValue(parameter,"polarisation") < 0.0) | (vectorGetValue(parameter,"polarisation") > pi)
      | (vectorGetValue(parameter,"phase") < 0.0) | (vectorGetValue(parameter,"phase") > 2.0*pi)
      | (vectorGetValue(parameter,"inclination") < 0.0) | (vectorGetValue(parameter,"inclination") > pi)){
    prior = -HUGE_VAL;
  }
  else {
    /* distance prior proportional to dl squared: */
    prior += 2.0 * vectorGetValue(parameter,"logdistance");
    prior += vectorIsElement(&MF->fixed,"latitude")    ? 0.0 : log(cos(vectorGetValue(parameter,"latitude")));
    prior += vectorIsElement(&MF->fixed,"inclination") ? 0.0 : log(sin(vectorGetValue(parameter,"inclination")));
    /* account for transformation of (mass1,mass2) --> (mc,eta): */
    prior += logJacobianMcEta(vectorGetValue(parameter,"chirpmass"), vectorGetValue(parameter,"massratio"));
    /* account for transformation  dL --> log(dL) : */
    prior += vectorGetValue(parameter,"logdistance");
    /* include the `detection probability' dependent on the signal amplitude: */
    ampli = (0.5*log(vectorGetValue(parameter,"massratio")) 
             + (5.0/6.0)*log(mc2mt(vectorGetValue(parameter,"chirpmass"),vectorGetValue(parameter,"massratio"))) 
             - vectorGetValue(parameter,"logdistance")
             + 0.5*log(1.0+cosiota2*(6.0+cosiota2)));
    prior -= log(1.0+exp((ampli-a)/b)); 
  }
  return prior;
}


double logpriorSineGaussianBurst(McmcFramework *MF, vector *parameter)
/* log - prior density for sine-Gaussian burst (8 parameters).      */
/* (not normalised!)                                                */
{
  double amp, sig;
  double Eamp, Esig;
  double prior = 0.0;
  /* check whether parameter values fall within domain: */
  if ((vectorGetValue(parameter,"frequency")<vectorGetValue(&MF->priorparameters,"freqLower")) 
      | (vectorGetValue(parameter,"frequency")>vectorGetValue(&MF->priorparameters,"freqUpper"))
      | (vectorGetValue(parameter,"time")<vectorGetValue(&MF->priorparameters,"timeLower")) | (vectorGetValue(parameter,"time")>vectorGetValue(&MF->priorparameters,"timeUpper"))
      | (vectorGetValue(parameter,"longitude") < -pi) | (vectorGetValue(parameter,"longitude") > pi)
      | (fabs(vectorGetValue(parameter,"latitude")) > pi/2.0)
      | (vectorGetValue(parameter,"polarisation")<0.0) | (vectorGetValue(parameter,"polarisation")>pi)
      | (vectorGetValue(parameter,"phase")<0.0) | (vectorGetValue(parameter,"phase")>2.0*pi)){
    prior = -HUGE_VAL;
  }
  else {
    amp = exp(vectorGetValue(parameter,"logamplitude"));
    sig = exp(vectorGetValue(parameter,"logsigma"));
    Eamp = vectorGetValue(&MF->priorparameters,"ampliExpect");
    Esig = vectorGetValue(&MF->priorparameters,"sigmaExpect");
    /* exponential prior distribution for amplitude:      */
    prior += -log(Eamp) - (1.0/(Eamp))*amp;
    /* exponential prior distribution for sigma:          */
    prior += -log(Esig) - (1.0/(Esig))*sig;
    /* prior propto sin(altitude) for altitude:           */
    prior += vectorIsElement(&MF->fixed,"latitude") ? 0.0 : log(cos(vectorGetValue(parameter,"latitude")));
    /* account for transformation  ampli --> log(ampli) : */
    prior += vectorGetValue(parameter,"logamplitude");
    /* account for transformation  sigma --> log(sigma) : */
    prior += vectorGetValue(parameter,"logsigma");
  }
  return prior;
}


double logguessdensity(McmcFramework *MF, vector *parameter)
/* 'wrapper' function for "guessing" density.               */
/* Checks the "MF->parameterset" character string and calls */
/* the corresponding density function.                      */
{
  double result = -HUGE_VAL;
  if (MF->parameterset == InspiralNoSpin) {
    if (MF->guessparameters == NULL)
      result = logpriorInspiralNospin(MF, parameter);
    else 
      result = logguessInspiralNospin(MF, parameter);
  }
  else if (MF->parameterset == BurstSineGaussian)
    result = logpriorSineGaussianBurst(MF, parameter);
  else
    printf(" : ERROR: referring to undefined parameterset in 'logguessdensity()'!\n");
  return(result);
}


double logguessInspiralNospin(McmcFramework *MF, vector *parameter)
/* log - "guessing" density for inspiral signal w/ no spin (9 parameters). */
/* (not normalised!)                                                       */
{
  double logdensity;
  /* check whether parameter values fall within domain: */
  if ((vectorGetValue(parameter,"massratio")>0.25) | (vectorGetValue(parameter,"massratio")<=0.0)
      | (mc2mass1(vectorGetValue(parameter,"chirpmass"), vectorGetValue(parameter,"massratio")) < vectorGetValue(&MF->priorparameters,"massLower"))
      | (mc2mass2(vectorGetValue(parameter,"chirpmass"), vectorGetValue(parameter,"massratio")) > vectorGetValue(&MF->priorparameters,"massUpper"))
      | (vectorGetValue(parameter,"time") < vectorGetValue(&MF->priorparameters,"timeLower")) 
      | (vectorGetValue(parameter,"time") > vectorGetValue(&MF->priorparameters,"timeUpper"))
      | (vectorGetValue(parameter,"longitude") < -pi) | (vectorGetValue(parameter,"longitude") > pi)
      | (fabs(vectorGetValue(parameter,"latitude")) > pi/2.0)
      | (vectorGetValue(parameter,"polarisation") < 0.0) | (vectorGetValue(parameter,"polarisation") > pi)
      | (vectorGetValue(parameter,"phase") < 0.0) | (vectorGetValue(parameter,"phase") > 2.0*pi)
      | (vectorGetValue(parameter,"inclination") < 0.0) | (vectorGetValue(parameter,"inclination") > pi)){
    logdensity = -HUGE_VAL;
  }
  else {
    /* log-normal distn. for mc: */
    logdensity = -log(vectorGetValue(parameter,"chirpmass")) 
      - pow(log(vectorGetValue(parameter,"chirpmass"))-log(MF->guessparameters[0]),2.0) / (2.0*0.01*0.01);
    /* normal for eta (including "bouncing off" the upper bound): */
    logdensity += log(exp(-pow(vectorGetValue(parameter,"massratio")-MF->guessparameters[1],2.0) / (2.0*0.05*0.05))
      +exp(-pow((0.25+(0.25-vectorGetValue(parameter,"massratio")))-MF->guessparameters[1],2.0) / (2.0*0.05*0.05)));
    /* normal for tc: */
    logdensity += -pow(vectorGetValue(parameter,"time")-MF->guessparameters[2],2.0) / (2.0*0.005*0.005);
    /* normal for log-distance: */
    logdensity += -pow(vectorGetValue(parameter,"logdistance")-MF->guessparameters[3],2.0) / (2.0*0.333*0.333);
    logdensity += vectorIsElement(&MF->fixed,"latitude")    ? 0.0 : log(cos(vectorGetValue(parameter,"latitude")));
    logdensity += vectorIsElement(&MF->fixed,"inclination") ? 0.0 : log(sin(vectorGetValue(parameter,"inclination")));
  }
  return logdensity;
}


void priordraw(McmcFramework *MF, vector *parameter)
/* 'wrapper' function for prior sampling.                   */
/* Checks the "MF->parameterset" character string and calls */
/* the corresponding prior sampling function.               */
{
  int i;
  if (MF->parameterset == InspiralNoSpin)
    priordrawInspiralNospin(MF, parameter);
  else if (MF->parameterset == BurstSineGaussian)
    priordrawBurstSineGaussian(MF, parameter);
  else
    printf(" : ERROR: referring to undefined parameterset in 'priordraw()'!\n");
  /* make sure fixed values remain unchanged: */
  for (i=0; i<parameter->dimension; ++i)
    if (vectorIsElement(&MF->fixed, parameter->name[i])) 
      vectorSetValue(parameter, parameter->name[i], vectorGetValue(&MF->fixed,parameter->name[i]));
}


void priordrawInspiralNospin(McmcFramework *MF, vector *parameter)
/* generates a random draw from the prior distribution. */
/* ...using a kind of rejection sampling: sampling from */
/* the `occurence' density is straightforward,          */
/* `undetected' samples are rejected...                 */
{
  double optiampli = log(0.5) + (5.0/6.0)*log(4.0) + 0.5*log(8.0);
  double ampli90 = optiampli - log(vectorGetValue(&MF->priorparameters,"dist90"));
  double ampli10 = optiampli - log(vectorGetValue(&MF->priorparameters,"dist10"));
  double a = (ampli10+ampli90)/2.0;
  double b = (ampli10-ampli90)/(-2.0*log(0.1/0.9));
  double ampli;
  double m1, m2, eta, logdist, cosiota, cosiota2;
  int detect = 0;
  double detectionprob;
  double maxdist = (vectorGetValue(&MF->priorparameters,"dist10")+2*(vectorGetValue(&MF->priorparameters,"dist10")-vectorGetValue(&MF->priorparameters,"dist90"))) 
                   * pow(vectorGetValue(&MF->priorparameters,"massUpper")/2.0, 5.0/6.0);
  while (!detect){
    m1            = gsl_ran_flat(GSLrandom, vectorGetValue(&MF->priorparameters,"massLower"), 
                           vectorGetValue(&MF->priorparameters,"massUpper"));
    m2            = gsl_ran_flat(GSLrandom, vectorGetValue(&MF->priorparameters,"massLower"), 
                           vectorGetValue(&MF->priorparameters,"massUpper"));
    eta           = (m1/(m1+m2)) * (m2/(m1+m2));
    logdist       = (1.0/3.0)*log(gsl_rng_uniform(GSLrandom)) + log(maxdist);
    cosiota       = gsl_ran_flat(GSLrandom, -1.0, 1.0);
    cosiota2      = cosiota * cosiota;
    ampli         = 0.5*log(eta) + (5.0/6.0)*log(m1+m2)
                    + 0.5 * log(1.0 + 6.0*cosiota2 + cosiota2*cosiota2)
                    - logdist;
    detectionprob = 1.0/(1.0+exp((ampli-a)/b));
    detect        = gsl_rng_uniform(GSLrandom) < detectionprob;
  }
  vectorSetValue(parameter, "chirpmass",    pow(m1*m2,0.6)/pow(m1+m2,0.2));
  vectorSetValue(parameter, "massratio",    eta);
  vectorSetValue(parameter, "logdistance",  logdist);
  vectorSetValue(parameter, "inclination",  acos(cosiota));
  vectorSetValue(parameter, "time",         gsl_ran_flat(GSLrandom, vectorGetValue(&MF->priorparameters,"timeLower"), 
                                                   vectorGetValue(&MF->priorparameters,"timeLower")));
  vectorSetValue(parameter, "longitude",    gsl_ran_flat(GSLrandom, -pi, pi));
  vectorSetValue(parameter, "latitude",     asin(gsl_ran_flat(GSLrandom, -1.0, 1.0)));
  vectorSetValue(parameter, "polarisation", gsl_ran_flat(GSLrandom, 0.0, pi));
  vectorSetValue(parameter, "phase",        gsl_ran_flat(GSLrandom, 0.0, 2.0*pi));
}


void priordrawBurstSineGaussian(McmcFramework *MF, vector *parameter)
/* generates a random draw from the prior distribution. */
{
  vectorSetValue(parameter, "frequency",    gsl_ran_flat(GSLrandom, vectorGetValue(&MF->priorparameters,"freqLower"), 
                                                   vectorGetValue(&MF->priorparameters,"freqUpper")));
  vectorSetValue(parameter, "logamplitude", log(gsl_ran_exponential(GSLrandom,vectorGetValue(&MF->priorparameters,"ampliExpect"))));
  vectorSetValue(parameter, "logsigma",     log(gsl_ran_exponential(GSLrandom,vectorGetValue(&MF->priorparameters,"sigmaExpect"))));
  vectorSetValue(parameter, "time",         gsl_ran_flat(GSLrandom, vectorGetValue(&MF->priorparameters,"timeLower"), 
                                                   vectorGetValue(&MF->priorparameters,"timeUpper")));
  vectorSetValue(parameter, "longitude",    gsl_ran_flat(GSLrandom, -pi, pi));
  vectorSetValue(parameter, "latitude",     asin(gsl_ran_flat(GSLrandom, -1.0,1.0)));
  vectorSetValue(parameter, "polarisation", gsl_ran_flat(GSLrandom, 0.0, pi));
  vectorSetValue(parameter, "phase",        gsl_ran_flat(GSLrandom, 0.0, 2.0*pi));
}


void guess(McmcFramework *MF, vector *parameter)
/* 'wrapper' function for sampling "guesses".               */
/* Checks the "MF->parameterset" character string and calls */
/* the corresponding "guessing" function.                   */
{
  int i;
  if (MF->parameterset == InspiralNoSpin){
    if (MF->guessparameters == NULL)
      priordrawInspiralNospin(MF, parameter);
    else
      guessInspiralNospin(MF, parameter);
  }
  else if (MF->parameterset == BurstSineGaussian)
    priordrawBurstSineGaussian(MF, parameter);
  else
    printf(" : ERROR: referring to undefined parameterset in 'guess()'!\n");
  /* make sure fixed values remain unchanged: */
  for (i=0; i<parameter->dimension; ++i)
    if (vectorIsElement(&MF->fixed, parameter->name[i])) 
      vectorSetValue(parameter, parameter->name[i], vectorGetValue(&MF->fixed,parameter->name[i]));
}

void guessInspiralNospin(McmcFramework *MF, vector *parameter)
/* Generates a random draw from the "guessing" distribution.  */
/* Note that changes in distribution parameters here need     */
/* to be accounted for in "logguessInspiralNospin()" as well! */
{
  double mc, eta, tc, logdist, m1, m2;
  int inRange = 0;
  while (!inRange){
    /* (note the -by now- hardcoded variances here) */
    mc      = exp(log(MF->guessparameters[0])+gsl_ran_ugaussian(GSLrandom)*0.01); /* mc  : 1%   */
    eta     = MF->guessparameters[1] + gsl_ran_ugaussian(GSLrandom)*0.05;         /* eta : 0.05 */
    tc      = MF->guessparameters[2] + gsl_ran_ugaussian(GSLrandom)*0.005;        /* tc  : 5 ms */
    logdist = log(MF->guessparameters[3]) + gsl_ran_ugaussian(GSLrandom)*0.333;   /* dist: 33%  */
    if (eta>0.25) eta = 0.25 - (eta-0.25);
    if (eta>0.0) {
      m1 = mc2mass1(mc, eta);
      m2 = mc2mass2(mc, eta);
      inRange = ((m1>vectorGetValue(&MF->priorparameters,"massLower")) 
                 && (m2<vectorGetValue(&MF->priorparameters,"massUpper")) 
                 && (tc>vectorGetValue(&MF->priorparameters,"timeLower"))
                 && (tc<vectorGetValue(&MF->priorparameters,"timeUpper")));
    }
  }
  vectorSetValue(parameter, "chirpmass",    mc);
  vectorSetValue(parameter, "massratio",    eta);
  vectorSetValue(parameter, "logdistance",  logdist);
  vectorSetValue(parameter, "inclination",  acos(gsl_ran_flat(GSLrandom, -1.0, 1.0)));
  vectorSetValue(parameter, "time",         tc);
  vectorSetValue(parameter, "longitude",    gsl_ran_flat(GSLrandom, -pi, pi));
  vectorSetValue(parameter, "latitude",     asin(gsl_ran_flat(GSLrandom, -1.0, 1.0)));
  vectorSetValue(parameter, "polarisation", gsl_ran_flat(GSLrandom, 0.0, pi));
  vectorSetValue(parameter, "phase",        gsl_ran_flat(GSLrandom, 0.0, 2.0*pi));
}


void propose(McmcFramework *MF, DataFramework *DF, int coherentN, 
             vector *parameter, double *logMHcoef)
/* 'wrapper' function for Metropolis-proposals.                  */
/* Checks the "MF->parameterset" character string,               */
/* calls the corresponding proposal function,                    */
/* and, based on the current MCMC state supplied in "parameter", */
/* returns a proposal, also in "parameter".                      */
/* "logMHcoef" is the coefficient making the difference between  */
/* the Metropolis- and Metropolis-Hastings algorithms:           */
/* the (log-) ratio of proposal densities                        */
/* for "backward" and "forward" proposals                        */
/*  J(theta^t-1 | theta^star)  /  J(theta^star | theta^t-1)      */
{
  int i;
  if (MF->parameterset == InspiralNoSpin)
    proposeInspiralNospin(MF, DF, coherentN, parameter, logMHcoef);
  else if (MF->parameterset == BurstSineGaussian)
    proposeBurstSineGaussian(MF, DF, coherentN, parameter, logMHcoef);
  else
    printf(" : ERROR: referring to undefined parameterset in 'propose()'!\n");
  /* double check fixed values remain unchanged: */
  for (i=0; i<parameter->dimension; ++i)
    if (vectorIsElement(&MF->fixed, parameter->name[i])) 
      vectorSetValue(parameter, parameter->name[i], vectorGetValue(&MF->fixed,parameter->name[i]));
}


void proposeInspiralNospin(McmcFramework *MF, DataFramework *DF, int coherentN, 
                           vector *parameter, double *logMHcoef)
/* Generate proposal for the MF->parameterset == "InspiralNoSpin" case,  */
/* i.e. when using the "20SP", "25SP", "2025" or "2535" waveform model.  */
{
  double *jump;
  double scalar;
  double dummy;
  int i, d;
  double skylocAngle, lati, longi;
  double posvec[3], pivotvec[3];

  /* at first, "MH-coefficient" (ratio of prop. densities) is one: */
  *logMHcoef = 0.0;
  /* (by now also remains at unity)                                */

  /* figure out actual parameter space dimension:                       */
  d = MF->startvalue.dimension - MF->fixed.dimension;
  /* adjust proposal scale accordingly:                                 */
  scalar = (d<=5) ? 1.0 : 2.4/sqrt(d);
  /* (see Gelman & al.: Bayesian Data Analysis, p.334)                  */
  scalar /= 5.0;
  /* (this seems to work better here than Gelman & al's recommendation) */

  /* generate the Gaussian / t-distributed proposal */
  /* using the random library (randlib) function:   */
  jump = (double*) malloc(sizeof(double) * MF->startvalue.dimension);
  RandMVNorm(MF, jump);
  if (MF->studentDF > 0.0)
    /* create Student-t out of normal                       */
    /* by multiplying with another (Chi^2-) Random Variable */
    for (i=0; i<MF->startvalue.dimension; ++i)
      jump[i] *= scalar * sqrt(MF->studentDF/gsl_ran_gamma(GSLrandom, MF->studentDF/2.0, 2.0));
    /* chisquared(nu) == gamma(nu/2,2)  (since GSL's chisq function seems buggy) */
  else /* Normal (Gaussian) proposal:                       */
    for (i=0; i<MF->startvalue.dimension; ++i)
      jump[i] *= scalar;

  /* invert every other inclination proposal: */
  if (gsl_rng_uniform(GSLrandom) < 0.5)
    jump[parindex(MF,"inclination")] *= -1.0;

  /* inflate some of proposals (by x4 or x4x4): */
  if (gsl_rng_uniform(GSLrandom) < 0.06)
    jump[parindex(MF,"chirpmass")] *= (gsl_rng_uniform(GSLrandom)<0.33) ? 4.0 : ((gsl_rng_uniform(GSLrandom)<0.5) ? 8.0 : 64.0);
  if (gsl_rng_uniform(GSLrandom) < 0.05)
    jump[parindex(MF,"massratio")] *= 4.0;

  /* Sky location proposals are treated specially (see below). */
  /* Store latitude jump, ignore longitude jump and do nothing yet. */
  skylocAngle = jump[parindex(MF,"latitude")];
  jump[parindex(MF,"latitude")] = 0.0;
  jump[parindex(MF,"longitude")] = 0.0;

  /* add 'jump' to current parameter value:                     */
  /* (i.e., proposal is -so far- centered around current value) */
  for (i=0; i<parameter->dimension; ++i)
    parameter->value[i] += jump[parindex(MF,parameter->name[i])];
  free(jump);

  /* manipulate sky location;                         */
  /* determine line-of-sight vector:                  */
  coord2vec(vectorGetValue(parameter,"latitude"), vectorGetValue(parameter,"longitude"), posvec);
  /* determine an (arbitrary) orthogonal unit vector: */
  pivotvec[0] = 0.0;
  pivotvec[1] = 1.0;
  pivotvec[2] =  - posvec[1] / posvec[2];
  normalise(pivotvec);
  /* rotate pivot into a random direction:  */
  rotate(pivotvec, gsl_ran_flat(GSLrandom, 0.0, 2.0*pi), posvec);
  /*  10% proposals in particular direction */
  /*  (-- iff 2-detector-network)           */
  if (gsl_rng_uniform(GSLrandom) < 0.10) { 
    /* check how many (different) detectors: */
    if (numberIfoSites(DF,coherentN)==2) {
      /* ifo #0 is the "1st" interferometer, now find second: */
      i = 1;
      while ((DF[i].ifo->latitude == DF[0].ifo->latitude) 
              && (DF[i].ifo->longitude == DF[0].ifo->longitude)) ++i;
      if (i >= coherentN) printf(" : WARNING: something's definitely fishy here!!\n");
      /* set "pivotvec" to be the connecting line between the two interferometers: */
      pivotvec[0] = DF[0].ifo->positionVector[0] - DF[i].ifo->positionVector[0];
      pivotvec[1] = DF[0].ifo->positionVector[1] - DF[i].ifo->positionVector[1];
      pivotvec[2] = DF[0].ifo->positionVector[2] - DF[i].ifo->positionVector[2];
      normalise(pivotvec);
    }
  }
  /* rotate position to get new position:  */
  rotate(posvec, skylocAngle, pivotvec);
  /* compute new (geographic) coordinates: */
  vec2coord(posvec, &lati, &longi);
  vectorSetValue(parameter, "latitude", lati);
  vectorSetValue(parameter, "longitude", longi);

  /* some uniform proposals: */
  if (gsl_rng_uniform(GSLrandom) < 0.10)  /* uniform mc proposal:    */
    vectorSetValue(parameter, "chirpmass",
                   gsl_ran_flat(GSLrandom, pow(vectorGetValue(&MF->priorparameters,"massLower"),1.2)/pow(2.0*vectorGetValue(&MF->priorparameters,"massLower"),0.2),
                          pow(vectorGetValue(&MF->priorparameters,"massUpper"),1.2)/pow(2.0*vectorGetValue(&MF->priorparameters,"massUpper"),0.2)));
  if (gsl_rng_uniform(GSLrandom) < 0.05)  /* uniform eta proposal:   */
    vectorSetValue(parameter, "massratio",
                   gsl_ran_flat(GSLrandom, (vectorGetValue(&MF->priorparameters,"massLower")*vectorGetValue(&MF->priorparameters,"massUpper"))/pow(vectorGetValue(&MF->priorparameters,"massLower")+vectorGetValue(&MF->priorparameters,"massUpper"),2.0),0.25));
  if (gsl_rng_uniform(GSLrandom) < 0.05)  /* uniform phase proposal: */
    vectorSetValue(parameter,"phase", gsl_ran_flat(GSLrandom, 0.0, 2.0*pi));
  if (gsl_rng_uniform(GSLrandom) < 0.05)  /* uniform t_c proposal:   */
    vectorSetValue(parameter,"time",
                   gsl_ran_flat(GSLrandom, vectorGetValue(&MF->priorparameters,"timeLower"), vectorGetValue(&MF->priorparameters,"timeUpper")));


  /* some proposals to related parameter space regions: */
  if (gsl_rng_uniform(GSLrandom)<0.05)   /* `pi-iota' - proposal: */
    vectorSetValue(parameter, "inclination", pi - vectorGetValue(parameter,"inclination"));
  if (gsl_rng_uniform(GSLrandom)<0.05)   /* `phi+pi' - proposal:  */
    vectorSetValue(parameter, "phase", vectorGetValue(parameter,"phase") + pi);
  if (gsl_rng_uniform(GSLrandom)<0.05)   /* `psi+pi' - proposal:  */
    vectorSetValue(parameter, "polarisation", vectorGetValue(parameter,"polarisation") + pi/2.0);


  /*--- "wrap around" out-of-domain proposals: ---*/
  dummy = vectorGetValue(parameter, "phase");
  if (dummy > 2.0*pi) {
    while (dummy > 2.0*pi) dummy -= 2.0*pi;
    vectorSetValue(parameter, "phase", dummy);
  }
  else if (dummy < 0.0) {
    while (dummy < 0.0) dummy += 2.0*pi;
    vectorSetValue(parameter, "phase", dummy);
  }

  dummy = vectorGetValue(parameter, "longitude");
  if (dummy > pi) {
    while (dummy > pi) dummy -= 2.0*pi;
    vectorSetValue(parameter, "longitude", dummy);
  }
  else if (dummy < -pi) {
    while (dummy < -pi) dummy += 2.0*pi;
    vectorSetValue(parameter, "longitude", dummy);
  }

  dummy = vectorGetValue(parameter, "polarisation");
  if (dummy > pi) {
    while (dummy > pi) dummy -= pi;
    vectorSetValue(parameter, "polarisation", dummy);
  }
  else if (dummy < 0.0) {
    while (dummy < 0.0) dummy += pi;
    vectorSetValue(parameter, "polarisation", dummy);
  }

  dummy = vectorGetValue(parameter, "latitude");
  if (dummy > pi/2.0){
    while (dummy > 1.5*pi) dummy -= 2.0*pi;
    if (dummy > pi/2.0) dummy -= 2.0*(dummy-pi/2.0);
    vectorSetValue(parameter, "latitude", dummy);
  }
  else if (dummy < -pi/2.0) {
    while (dummy < -1.5*pi) dummy += 2.0*pi;
    if (dummy < -pi/2.0) dummy += 2.0*(-pi/2.0-dummy);
    vectorSetValue(parameter, "latitude", dummy);
  }

  dummy = vectorGetValue(parameter, "inclination");
  if (dummy > pi){
    while (dummy > 2.0*pi) dummy -= 2.0*pi;
    if (dummy > pi) dummy -= 2.0*(dummy-pi);
    vectorSetValue(parameter, "inclination", dummy);
  }
  else if (dummy < 0.0) {
    while (dummy < -1.0*pi) dummy += 2.0*pi;
    if (dummy < 0.0) dummy += -2.0*dummy;
    vectorSetValue(parameter, "inclination", dummy);
  }

  if (vectorGetValue(parameter,"massratio") > 0.25) 
    vectorSetValue(parameter,"massratio", 0.5-vectorGetValue(parameter,"massratio"));
  /*--- (end wrapping) ---*/
}


void proposeBurstSineGaussian(McmcFramework *MF, DataFramework *DF, int coherentN,
                              vector *parameter, double *logMHcoef)
/* Generate proposal for the  MF->parameterset == BurstSineGaussian  case. */
{
  double *jump;
  double scalar;
  double dummy;
  int i, d;
  double skylocAngle, lati, longi;
  double posvec[3], pivotvec[3];

  /* figure out actual parameter space dimension: */
  d = MF->startvalue.dimension - MF->fixed.dimension;
  /* adjust proposal scale accordingly: */
  scalar = (d<=5) ? 1.0 : 2.4/sqrt(d);
  /* (see Gelman & al.: Bayesian Data Analysis, p.334) */
  scalar /= 10.0;

  /* generate the Gaussian / t-distributed proposal */
  /* using the random library (randlib) function:   */
  jump = (double*) malloc(sizeof(double) * MF->startvalue.dimension);
  RandMVNorm(MF, jump);
  if (MF->studentDF > 0.0)
    /* create Student-t out of normal                       */
    /* by multiplying with another (Chi^2-) Random Variable */
    for (i=0; i<MF->startvalue.dimension; ++i)
      jump[i] *= scalar * sqrt(MF->studentDF/gsl_ran_gamma(GSLrandom, MF->studentDF/2.0, 2.0));
    /* chisquared(nu) == gamma(nu/2,2)  (since GSL's chisq function seems buggy) */
  else /* Normal (Gaussian) proposal: */
    for (i=0; i<MF->startvalue.dimension; ++i)
      jump[i] *= scalar;

  /* inflate some of proposals: */
  if (gsl_rng_uniform(GSLrandom)<0.1)
    jump[parindex(MF,"frequency")] *= (gsl_rng_uniform(GSLrandom)<0.5) ? 10.0 : 100.0;

  /* Sky location proposals are treated specially (see below).      */
  /* Store latitude jump, ignore longitude jump and do nothing yet. */
  skylocAngle = jump[parindex(MF,"latitude")];
  jump[parindex(MF,"latitude")] = 0.0;
  jump[parindex(MF,"longitude")] = 0.0;

  /* add 'jump' to current parameter value:                     */
  /* (i.e., proposal is -so far- centered around current value) */
  for (i=0; i<parameter->dimension; ++i)
    parameter->value[i] += jump[parindex(MF,parameter->name[i])];
  free(jump);

  /* manipulate sky location;                         */
  /* determine line-of-sight vector:                  */
  coord2vec(vectorGetValue(parameter,"latitude"), vectorGetValue(parameter,"longitude"), posvec);
  /* determine an (arbitrary) orthogonal unit vector: */
  pivotvec[0] = 0.0;
  pivotvec[1] = 1.0;
  pivotvec[2] =  - posvec[1] / posvec[2];
  normalise(pivotvec);
  /* rotate pivot into a random direction:  */
  rotate(pivotvec, gsl_ran_flat(GSLrandom, 0.0, 2.0*pi), posvec);
  /*  10% proposals in particular direction */
  /*  (-- iff 2-detector-network)           */
  if (gsl_rng_uniform(GSLrandom) < 0.10) { 
    /* check how many (different) detectors: */
    if (numberIfoSites(DF,coherentN)==2) {
      /* ifo #0 is the "1st" interferometer, now find second: */
      i = 1;
      while ((DF[i].ifo->latitude == DF[0].ifo->latitude) 
              && (DF[i].ifo->longitude == DF[0].ifo->longitude)) ++i;
      if (i >= coherentN) printf(" : WARNING: something's definitely fishy here!!\n");
      /* set "pivotvec" to be the connecting line between the two interferometers: */
      pivotvec[0] = DF[0].ifo->positionVector[0] - DF[i].ifo->positionVector[0];
      pivotvec[1] = DF[0].ifo->positionVector[1] - DF[i].ifo->positionVector[1];
      pivotvec[2] = DF[0].ifo->positionVector[2] - DF[i].ifo->positionVector[2];
      normalise(pivotvec);
    }
  }
  /* rotate position to get new position: */
  rotate(posvec, skylocAngle, pivotvec);
  /* compute new (geographic) coordinates: */
  vec2coord(posvec, &lati, &longi);
  vectorSetValue(parameter, "latitude", lati);
  vectorSetValue(parameter, "longitude", longi);

  /* some uniform proposals: */
  if (gsl_rng_uniform(GSLrandom) < 0.05) /* uniform frequency proposal: */
    vectorSetValue(parameter,"frequency",
                   gsl_ran_flat(GSLrandom, vectorGetValue(&MF->priorparameters,"freqLower"), 
                          vectorGetValue(&MF->priorparameters,"freqUpper")));
  if (gsl_rng_uniform(GSLrandom)<0.05)
    vectorSetValue(parameter,"phase", gsl_ran_flat(GSLrandom, 0.0, 2.0*pi));
  if (gsl_rng_uniform(GSLrandom)<0.05)
    vectorSetValue(parameter, "time",
                   gsl_ran_flat(GSLrandom, vectorGetValue(&MF->priorparameters,"timeLower"), 
                          vectorGetValue(&MF->priorparameters,"timeUpper")));

  /* time/phase: */
  if (gsl_rng_uniform(GSLrandom) < 0.05) {
    vectorSetValue(parameter, "time",
                   vectorGetValue(parameter,"time")
                   + gsl_ran_ugaussian(GSLrandom)*0.5*exp(vectorGetValue(parameter,"logsigma")));
    vectorSetValue(parameter,"phase", gsl_ran_flat(GSLrandom, 0.0, 2.0*pi));
  }
  /* uniform extrinsic parameters: */
  if (gsl_rng_uniform(GSLrandom) < 0.2) {
    vectorSetValue(parameter, "time",
                   gsl_ran_flat(GSLrandom, vectorGetValue(&MF->priorparameters,"timeLower"), 
                   vectorGetValue(&MF->priorparameters,"timeUpper")));
    vectorSetValue(parameter,"phase", gsl_ran_flat(GSLrandom, 0.0, 2.0*pi));    
    vectorSetValue(parameter,"polarisation", gsl_ran_flat(GSLrandom, 0.0, pi));
    vectorSetValue(parameter,"latitude", gsl_ran_flat(GSLrandom, -pi/2.0, pi/2.0));
    vectorSetValue(parameter,"longitude", gsl_ran_flat(GSLrandom, -pi, pi));
  }

  /* some proposals to related parameter space regions: */
  if (gsl_rng_uniform(GSLrandom)<0.05)   /* `phi+pi' - proposal: */
    vectorSetValue(parameter, "phase", vectorGetValue(parameter,"phase")+pi);
  if (gsl_rng_uniform(GSLrandom)<0.05)   /* `tc+1/f' - proposal */
    vectorSetValue(parameter, "time",
                   vectorGetValue(parameter,"time") + (gsl_rng_uniform(GSLrandom)<0.5 ? -1.0/vectorGetValue(parameter,"frequency") : 1.0/vectorGetValue(parameter,"frequency")));
  if (gsl_rng_uniform(GSLrandom)<0.05)   /* `psi+pi/2' - proposal: */
    vectorSetValue(parameter, "polarisation", vectorGetValue(parameter,"polarisation")+pi/2.0);

  if (gsl_rng_uniform(GSLrandom)<0.05){   /* `psi+pi/2' AND `phi+pi'- proposal: */
    vectorSetValue(parameter, "polarisation", vectorGetValue(parameter,"polarisation")+pi/2.0);
    vectorSetValue(parameter, "phase", vectorGetValue(parameter,"phase")+pi);
  }


  /*--- "wrap around" out-of-domain proposals: ---*/
  dummy = vectorGetValue(parameter, "phase");
  if (dummy > 2.0*pi) {
    while (dummy > 2.0*pi) dummy -= 2.0*pi;
    vectorSetValue(parameter, "phase", dummy);
  }
  else if (dummy < 0.0) {
    while (dummy < 0.0) dummy += 2.0*pi;
    vectorSetValue(parameter, "phase", dummy);
  }

  dummy = vectorGetValue(parameter, "longitude");
  if (dummy > pi) {
    while (dummy > pi) dummy -= 2.0*pi;
    vectorSetValue(parameter, "longitude", dummy);
  }
  else if (dummy < -pi) {
    while (dummy < -pi) dummy += 2.0*pi;
    vectorSetValue(parameter, "longitude", dummy);
  }

  dummy = vectorGetValue(parameter, "polarisation");
  if (dummy > pi) {
    while (dummy > pi) dummy -= pi;
    vectorSetValue(parameter, "polarisation", dummy);
  }
  else if (dummy < 0.0) {
    while (dummy < 0.0) dummy += pi;
    vectorSetValue(parameter, "polarisation", dummy);
  }

  dummy = vectorGetValue(parameter, "latitude");
  if (dummy > pi/2.0){
    while (dummy > 1.5*pi) dummy -= 2.0*pi;
    if (dummy > pi/2.0) dummy -= 2.0*(dummy-pi/2.0);
    vectorSetValue(parameter, "latitude", dummy);
  }
  else if (dummy < -pi/2.0) {
    while (dummy < -1.5*pi) dummy += 2.0*pi;
    if (dummy < -pi/2.0) dummy += 2.0*(-pi/2.0-dummy);
    vectorSetValue(parameter, "latitude", dummy);
  }
  /*--- (end wrapping) ---*/


  *logMHcoef = 0.0;
}


void importanceresample(DataFramework *DF, int coherentN, McmcFramework *MF,
                        vector *parameter, 
                        long samplesize, long subsamplesize)
/* generates a single draw from the posterior distribution */
/* using importance resampling                             */
/* (see e.g. Gelman et al: Bayesian Data Analysis).        */

/* NOTE: for generality,                                   */
/*  "parameter" is assumed to be                           */
/* a vector of 'subsamplesize' elements of type 'vector',  */
/*  though for now 'subsamplesize' will usually be ==1.    */
{
  struct listelement{
    vector             *par;
    double             logposterior;
    double             impratio;
    struct listelement *next;
  };
  struct listelement *list=NULL;
  struct listelement *new=NULL;
  struct listelement *search=NULL;
  int  listlength=0, recycle = 0;
  int  i, j;
  double maximp, impsum, random, probsum;
  time_t starttime, endtime;
  double seconds;
  int n = samplesize;
  int m = subsamplesize;
  double postdiff = ((double)MF->startvalue.dimension)/2.0 
                    + 3.0*sqrt(((double)MF->startvalue.dimension)/2.0);
  /* n : number of prior draws (internal)                         */
  /* m : number of importance from above draws (result)           */
  /* postdiff : max. posterior difference                         */

  if (verbose) {
    printf(" | importance resampling: %d out of %d,\n", m, n);
    printf(" | posterior difference threshold = %.2f\n", postdiff);
  }

  if (m>n) printf(" : ERROR: Need to set 'subsamplesize' << 'samplesize' in 'importancesample()'!\n");

  time(&starttime);
  for (i=0; i<n; ++i) {
    if (!recycle) {
      new = (struct listelement*) malloc(sizeof(struct listelement));
      new->par = (vector*) malloc(sizeof(vector));
      vectorInit(new->par);
      vectorSetup(new->par, MF->parameterset);
    }
    guess(MF, new->par);
    new->logposterior  = logprior(MF, new->par);
    new->logposterior += loglikelihood(DF, coherentN, MF->template, new->par);
    new->next = NULL;
    /* if list short enough or new element's posterior large enough, then include in list: */
    if ((listlength < m) || ((new->logposterior) >= ((list->logposterior)-postdiff))){
      /*printf(" :      -> insert\n");*/
      if (list == NULL) {list = new; ++listlength;} /* insert as very very first element */
      else{
        if ((new->logposterior) > (list->logposterior)){ /* new element is better than best by now */
          /* insert as first element: */
          new->next = list;
          list = new;
          ++listlength;
          /* search for elements to be thrown out: */
          search = list;
          j=1;
          while ((search->next != NULL) && ((search->next->logposterior >= ((list->logposterior)-postdiff))
                  || (j<m))){
            search = search->next;
            ++j;
          }
          /* now `search->next' and following ones are thrown out. */
          if (search->next != NULL){
            new = search->next;
            search->next = NULL;
            search = new->next;
            while (new != NULL){
	      /*printf(" :      -> pop (%f)\n",new->logposterior);*/
              vectorDispose(new->par);
              free(new->par);
              free(new);
              new = search;
              search = (search == NULL) ? NULL : search->next;
              --listlength;
            }
          }
        }
        else { /* no improvement, insert in list: */
          search = list;
          /* search place to insert: */
          while ((search->next != NULL) && (search->next->logposterior > new->logposterior))
            search = search->next;
          /* bend pointers: */
          new->next = search->next;
          search->next = new;
          ++listlength;
        }
      }
      recycle = 0;
    }
    else {
      recycle = 1;
      /*printf(" :      -> ignore\n");*/
    }
    if (((i+1)%5000==0) | (i==19) | (i==49) | (i==99) | (i==199) | (i==499) | (i==749) | (i==999) | (i==1999))
      if (verbose) printf(" | sample: %6d;  %5d in queue (top candidate: %.6f)\n", i+1, listlength, list->logposterior);
  }
  /*printf(" : drawing FINISHED\n");*/


  /* now replace loglikelihoods by log - importance ratios: */
  search = list;
  for (i=0; i<listlength; ++i){
    search->logposterior -= logguessdensity(MF, search->par);
    search = search->next;
  }

  /* draw remaining samples 1,...,(m-1) */
  for (j=0; j<m; ++j){
    /* determine maximum, shift values & compute actual ratios */
    maximp = -HUGE_VAL;
    search = list;
    for (i=1; i<=listlength; ++i){
      maximp = (search->logposterior > maximp) ? search->logposterior : maximp;  
      search = search->next;
    }
    impsum = 0.0;
    search = list;
    for (i=1; i<=listlength; ++i){
      impsum += search->impratio = exp(search->logposterior - maximp);
      search = search->next;
    }
    /* draw a parameter set: */
    random = gsl_ran_flat(GSLrandom, 0.0, impsum);
    search = list;
    probsum = list->impratio;
    while (probsum < random) {
      search = search->next;
      probsum += search->impratio;
    }
    /* copy drawn parameter set to result vector: */
    vectorCopy(search->par, &parameter[j]);
    /* exclude from further draws: */
    search->logposterior = -HUGE_VAL;
  }
  time(&endtime);
  seconds = difftime(endtime, starttime);
  if (verbose) printf(" | ...finished at %.3f samples per second.\n", ((double)n)/seconds);
  MF->secPerIteration = seconds/((double)n);
}


void logtofile(McmcFramework *MF, vector *parameter, 
               long iteration, long accepted, 
               double logprior, double loglikelihood, double logposterior)
/* "iteration" indicates the iteration number.                                  */
/* For  "iteration==LONG_MIN"  a new, empty file (header only) is created       */
/* (or possibly overwritten), otherwise a line is appended to an existing file. */
{
  FILE *logfile;
  int i;
  vector outvector; /* vector that is actually logged */

  /* determine what exactly is going to be logged:      */
  vectorInit(&outvector);
  if (MF->celestialCoords) {
    /* copy everything EXCEPT latitude / longitude,     */
    /* transform geographical to celestial coordinates: */
    for (i=0; i<parameter->dimension; ++i)
      if ((strcmp(parameter->name[i],"latitude")==0) && (!vectorIsElement(&MF->fixed, "latitude")))
        vectorAdd(&outvector, "declination", vectorGetValue(parameter, "latitude"));
      else if ((strcmp(parameter->name[i],"longitude")==0)  && (!vectorIsElement(&MF->fixed, "longitude")))
        vectorAdd(&outvector, "rightascension", 
                  rightAscension(vectorGetValue(parameter, "longitude"), MF->GMST));
      else
        vectorAdd(&outvector, parameter->name[i], parameter->value[i]);
  }  
  else vectorCopy(parameter, &outvector);

  if (iteration==LONG_MIN) {  /* set up new file: */
    if (verbose) printf(" | logging to file \"%s\"\n",MF->logfilename);
    logfile = fopen(MF->logfilename, "w");
    fprintf(logfile, "iteration accepted logprior loglikelihood logposterior");
    for (i=0; i<outvector.dimension; ++i)
      if ((!vectorIsElement(&MF->fixed, outvector.name[i]))
          || ((strcmp(outvector.name[i],"declination")==0) && (!vectorIsElement(&MF->fixed, "latitude")))
          || ((strcmp(outvector.name[i],"rightascension")==0) && (!vectorIsElement(&MF->fixed, "longitude"))))
        fprintf(logfile, " %s", outvector.name[i]);
    fprintf(logfile, "\n");
    fclose(logfile);
  }
  else { /* append to existing file: */
    logfile = fopen(MF->logfilename, "a");
    fprintf(logfile, "%d %d %.6f %.6f %.6f", iteration, accepted, 
            logprior, loglikelihood, logposterior);
    for (i=0; i<outvector.dimension; ++i)
      if (!vectorIsElement(&MF->fixed, outvector.name[i]))
        fprintf(logfile, " %.14e", outvector.value[i]);
    fprintf(logfile, "\n");
    fclose(logfile);
  }
  vectorDispose(&outvector);
}


void printtime()
/* prints time (& date) to screen */
{
  time_t tm;
  struct tm *ltime;
  time( &tm );
  ltime = localtime( &tm );
  ltime->tm_mon++;
  ltime->tm_year += 1900;
  printf("%02i.%02i.%04i  %02i:%02i:%02i\n", ltime->tm_mday, ltime->tm_mon,
          ltime->tm_year, ltime->tm_hour, ltime->tm_min, ltime->tm_sec);
}

void savePSD(DataFramework *DF, char *filename)
/* save the (1-sided) power spectral density (given in 'DF->powspec') to a file */
{
  int i;
  FILE *specfile;
  printf(" : writing (logarithmic, one-sided) noise power spectral density to file\n : '%s'...\n",
         filename);
  specfile = fopen(filename, "w");
  fprintf(specfile, "f logPSD\n");
  for (i=0; i<DF->FTSize; ++i)
    fprintf(specfile, "%f %e\n", ((double)i)*DF->FTDeltaF, DF->powspec[i]);
  fclose(specfile);
  printf(" : ...finished.\n");
}


void metropolishastings(McmcFramework *MF, DataFramework *DF, int coherentN)
/* The actual Metropolis (MCMC) sampler. */
{
  long i, acceptcount=0;
  vector state;             /* current state of MC                         */
  double lprior, llikeli;   /* corresponding prior & likelihood            */
  vector proposal;          /* new, proposed state of MC                   */
  double lpriorP, llikeliP; /* corresponding prior & likelihood            */
  double logMHcoef;         /* Metropolis-Hastings "asymmetry" coefficient */
  double logalpha;          /* log of acceptance probability               */
  int accept;
  time_t starttime, endtime;
  double seconds;

  vectorInit(&state);
  vectorSetup(&state, MF->parameterset);
  vectorInit(&proposal);
  vectorSetup(&proposal, MF->parameterset);

  /* create log file:                                                 */
  logtofile(MF, &state, LONG_MIN, 0, 0.0, 0.0, 0.0);

  /* first compute "null" likelihood, for no signal present:          */
  llikeli = loglikelihood(DF, coherentN, MF->template, NULL);
  /* log as "negative 1st" iteration, with all zero parameter values: */
  if (verbose) printf(" | 'null' likelihood         :  %.3f\n", llikeli);
  logtofile(MF, &state, -1, 0, 0.0, llikeli, 0.0);

  /* now proceed at starting parameter values:                        */
  vectorCopy(&MF->startvalue, &state);
  lprior  = logprior(MF, &state);
  if (lprior == -HUGE_VAL) {
    printf(" : WARNING: starting value has zero prior probability!\n");
    vectorPrint(&state);
  }
  llikeli = loglikelihood(DF, coherentN, MF->template, &state);
  if (verbose) printf(" | starting value likelihood :  %.3f\n", llikeli);
  logtofile(MF, &state, 0, 0, lprior, llikeli, lprior+llikeli);

  if (verbose){
    printf(" | starting Metropolis-sampler; signal template used: '#%d'\n | iterations: %d", 
           MF->template, MF->iterations);
    if (MF->secPerIteration > 0.0)
      printf(", estimated time: %.2f hrs", (MF->iterations*MF->secPerIteration)/3600.0);
    printf("\n");
  }
  if (verbose) {
    printf(" | started  : ");
    printtime();
  }

  /*-- start actual M-H sampler: --*/
  time(&starttime);
  for (i=1; i<=MF->iterations; ++i) {
    vectorCopy(&state, &proposal);
    propose(MF, DF, coherentN, &proposal, &logMHcoef);
    lpriorP  = logprior(MF, &proposal);
    if (lpriorP > -HUGE_VAL) { /* (only compute likelihood if prior wasn't already zero) */
      llikeliP = loglikelihood(DF, coherentN, MF->template, &proposal);
      logalpha = (lpriorP-lprior) + (llikeliP-llikeli) + logMHcoef;
      accept = ((logalpha>=0.0) || (log(gsl_rng_uniform(GSLrandom))<logalpha));
    }
    else {                     /* (don't bother computing likelihood, reject right away) */
      llikeliP = -HUGE_VAL;
      logalpha = -HUGE_VAL;
      accept   = 0;
    }
    if (accept){               /* update state                                           */
      vectorCopy(&proposal, &state);
      lprior  = lpriorP;
      llikeli = llikeliP;
      acceptcount += 1;
    }
    if (i % 100 == 0)          /* log every 100th iteration                              */
      logtofile(MF, &state, i, acceptcount, lprior, llikeli, lprior+llikeli);
    if (verbose && (fmod(log((double)i)/log(2.0),1.0) == 0.0)){
      time(&endtime);
      seconds = difftime(endtime, starttime);
      if (seconds > 0.0)
        printf(" | %.5f seconds per iteration (%.4f iter/sec).\n", 
               seconds/((double)i), ((double)i)/(seconds));
    }
  }
  time(&endtime);
  seconds = difftime(endtime, starttime);
  if (verbose) {
    printf(" | finished : ");
    printtime();
    printf(" |            (%.1f seconds)\n", seconds);
  }
  vectorDispose(&state);
  vectorDispose(&proposal);
}


void printDF(DataFramework *DF)
{
  int i;
  printf(" : ------------------------------------------------\n");
  printf(" : timeCenter                  %.3fs\n", DF->timeCenter);
  printf(" : timeBefore, -After          %.3fs, %.3f s\n", DF->timeBefore, DF->timeAfter);
  printf(" : minF-maxF                   %.1f-%.1f Hz\n", DF->minF, DF->maxF);
  printf(" : tukeypar                    %.3f\n", DF->tukeypar);
  printf(" : frameChannel                \"%s\"\n", DF->frameChannel);
  printf(" : datacachefile               \"%s\"\n", DF->datacachefile);
  printf(" : noisecachefile              \"%s\"\n", DF->noisecachefile);
  printf(" : datafilenames               \"%s\"\n", DF->datafilenames);
  printf(" : noisefilenames              \"%s\"\n", DF->noisefilenames);
  printf(" : PsdEstimateStart, -End, -N  %.3f, %.3f, %d\n", DF->PsdEstimateStart, DF->PsdEstimateEnd, DF->PsdEstimateN);
  printf(" : simulateData                %d\n", DF->simulateData);
  printf(" : simulatePsd                 %d\n", DF->simulatePsd);
  printf(" : dataStart                   %.3f\n", DF->dataStart);
  printf(" : dataDeltaT                  %.5f s\n", DF->dataDeltaT);
  printf(" : dataSize                    %d\n", DF->dataSize);
  printf(" : *data                       %d\n", DF->data);
  printf(" : *window                     %d\n", DF->window);
  printf(" : winss                       %.1f\n", DF->winss);
  printf(" : *dataFT                     %d\n", DF->dataFT);
  printf(" : FTsize                      %d\n", DF->FTSize);
  printf(" : FTDeltaF                    %.3f Hz\n", DF->FTDeltaF);
  printf(" : *powspec                    %d\n", DF->powspec);
  printf(" : noiseModel                  \"%s\"\n", DF->noiseModel);
  printf(" : minInd, maxInd              %d, %d\n", DF->minInd, DF->maxInd);
  printf(" : rawDataRange                %.3f s\n", DF->rawDataRange);
  printf(" : ifo->name                   \"%s\"\n", DF->ifo->name);
  if (DF->data != NULL) for (i=0; i<5; ++i) 
    printf(" :   data[%d] = %e\n", i, DF->data[i]);
  if (DF->powspec != NULL) for (i=0; i<5; ++i) 
    printf(" :   powspec[%d] = %e\n", i, DF->powspec[i]);
}


int numberIfoSites(DataFramework *DF, int coherentN)
{
  int i, j, found, count = 0;
  for (i=0; i<coherentN; ++i) {
    if (i==0) ++count;
    else {
      j = 0;
      found = 0;
      while ((!found) && (j<i)) {
        if ((DF[i].ifo->latitude != DF[j].ifo->latitude) 
            && (DF[i].ifo->longitude != DF[j].ifo->longitude)) found = 1;
	++j;
      }
      if (found) ++count;
    }
  }
  return count;
}

double GMST(double GPSsec)
/* Derives the `Greenwich Mean Sidereal Time' (in radians!) */
/* from GPS time (in seconds).                              */
/* (see K.R.Lang(1999): Astrophysical formulae, p.80 sqq.)  */
{
  /* double Julian_Jan1st2000midnight = 2451544.5; */
  /* double Julian_Jan1st2000noon     = 2451545.0; */
  double GPS_Jan1st2000midnight    = 630720013.0;
  double leapseconds = 32.0;  /* at Jan 1st 2000 */
  double seconds, days, centuries, secCurrentDay, result;
  if (GPSsec > (GPS_Jan1st2000midnight + 189388800.0)) 
    leapseconds += 1.0; /* one more leapsecond after 2005/'06 */
  if(GPSsec > 914803214.0) leapseconds += 1.0; /* Leap second after 2008/'09 */
  if (GPSsec < GPS_Jan1st2000midnight) {
    printf(" : WARNING: GMST's before 1.1.2000 may be inaccurate! \n");
    printf(" :          (requested: GMST(GPS=%.3fs))\n", GPSsec);
  }
  /* time since Jan 1st 2000 (0:00h) */
  seconds       = (GPSsec - GPS_Jan1st2000midnight) + (leapseconds - 32.0);
  days          = floor( seconds / (24.0*60.0*60.0) ) -0.5;
  secCurrentDay = fmod(seconds, 24.0*60.0*60.0);
  centuries     = days / 36525.0;
  result  = 24110.54841+(centuries*(8640184.812866+centuries*(0.093104+centuries*6.2e-6)));
  result += secCurrentDay * 1.002737909350795; /* (UTC day is 1.002 * MST day) */
  result  = fmod(result/(24.0*60.0*60.0),1.0);
  result *= 2.0*pi;
  return result;
}

double rightAscension(double longi, double gmst)
/* Derives right ascension (in radians!) from longitude given GMST (radians). */
/* Declination == latitude for equatorial coordinates.                        */
{
  double result = longi + gmst;
  while (result<0.0) result += 2.0*pi;
  while (result>2.0*pi) result -= 2.0*pi;
  /*result *= 24.0/(2.0*pi);*/
  return result;
}

double longitude(double rightascension, double gmst)
/* Derives longitude from right ascension (radians), given GMST (radians).    */
{
  double result = rightascension - gmst;
  while (result > pi)  result -= 2.0*pi;
  while (result < -pi) result += 2.0*pi;
  return result;
}


int main(int argc, char *argv[])
{
  DataFramework *DatFW;
  McmcFramework *McmcFW;
  interferometer *ifodata = NULL;
  int ifoN=0;    /* number of Ifos for which geographic location etc is available. */
  int coherentN; /* number of data sets (coherently) used.                         */
  int initOK;

  if (verbose) 
    printf(" +----[ followupMcmc ]---------------------------------------------------------\n");

  /* initialise GSL random number generator: */
  gsl_rng_env_setup();
  GSLrandom = gsl_rng_alloc(gsl_rng_mt19937);  /* use `Mersenne Twister' algorithm         */
  /* printf(" : Using GSL's \"%s\" random number generator.\n", gsl_rng_name (GSLrandom)); */
  /* (random seed is set later in `init()' function)                                       */

  if (argc==1)
    printhelpmessage();
  else {

    ifoInit(&ifodata, &ifoN);

    initOK = init(&DatFW, &McmcFW, &coherentN, argc, argv, ifodata, ifoN);

    /*
     *  savePSD(&DatFW[0], "/home/christian/temp/spec8.txt");
     *  initOK = 0;
     */

    /* printf(" : number of (different) ifo sites: %d\n", numberIfoSites(DatFW,coherentN)); */

    if (initOK) {
      /* for (i=0; i<coherentN; ++i) printDF(&DatFW[i]); */
      metropolishastings(McmcFW, DatFW, coherentN);
    }
    /* clearMF(McmcFW);          */
    /* clearDF(DatFW);           */
    /* clearIfo(ifodata, &ifoN); */
  }
  gsl_rng_free(GSLrandom);
  if (verbose) 
    printf(" +-----------------------------------------------------------------------------\n");

  return 0;
}
