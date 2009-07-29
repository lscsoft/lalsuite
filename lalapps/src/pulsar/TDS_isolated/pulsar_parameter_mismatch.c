#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALDatatypes.h>
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
#include <lal/LALRCSID.h>

#include <getopt.h>

#define SIXTH 0.166666666666666666666666666666666666666667L
#define TWENTYFOURTH 0.04166666666666666666666666666666666666666667L
#define MAXPARAMS 35

#define USAGE \
"Usage: %s [options]\n\n"\
" --help              display this message\n"\
" --detector          IFOs to use e.g. H1\n"\
" --pulsar            name of pulsar e.g. J0534+2200\n"\
" --par-file          pulsar parameter (.par) file (full path)\n"\
" --par-file2         pulsar parameter file to compare\n"\
" --covariance        pulsar covariance matrix file\n"\
" --start             start time GPS\n"\
" --timespan          time span to calculate over (seconds)\n"\
" --deltat            length of each step (seconds)\n"\
" --iterations        number of random points to check\n"\
" --output-dir        output directory\n"\
" --earth             Earth ephemeris\n"\
" --sun               Sun ephemeris\n"\
"\n"


RCSID("$Id$");

typedef struct tagInputParams{
  CHAR det[10];
  CHAR pulsar[15];
  CHAR parfile[256];
  CHAR parfile2[256];
  CHAR covfile[256];

  INT4 start;
  INT4 timespan;
  INT4 deltat;

  INT4 iterations;

  CHAR outputdir[256];
  CHAR earth[256];
  CHAR sun[256];
}InputParams;

typedef struct tagParamData{
  CHAR *name;  /* parameter name as given by the conventions in the .mat file
(see param variable in TEMPOs mxprt.f file */
  REAL8 val;   /* parameter value */
  REAL8 sigma; /* standard deviation on the parameter as read from the .par
file */
  INT4 matPos;  /* row position in the covariance matrix */
}ParamData;

void get_input_args(InputParams *inputParams, INT4 argc, CHAR *argv[]);

REAL8Array *CholeskyDecomp( REAL8Array *M, CHAR* uOrl );

REAL8Array *ReadCorrelationMatrix( CHAR *matrixFile, 
  BinaryPulsarParams params, ParamData *data );

REAL8Array *CreateCovarianceMatrix( ParamData *data, REAL8Array *corMat );

REAL8Array *XLALCheckPositiveDefinite( REAL8Array *matrix );

REAL8Array *XLALConvertToPositiveDefinite( REAL8Array *nonposdef );

ParamData *MultivariateNormalDeviates( REAL8Array *covmat, ParamData *means,
  RandomParams *randomParams );

void SetMCMCPulsarParams( BinaryPulsarParams *pulsarParams, ParamData *data );

REAL8Vector *get_phi( double start, double deltaT, int npoints,
 BinaryPulsarParams params, BarycenterInput bary, EphemerisData *edat );

int main(int argc, char *argv[]){
  static LALStatus status;

  InputParams inputs;

  FILE *fp=NULL;

  BinaryPulsarParams params, params2;

  BinaryPulsarParams deltaparams;

  int i=0, count=0, countBin=0, j=1, k=0, N=0;

  int npoints=0; /* number of 10 minutes stretches */

  EphemerisData *edat=NULL;
  LALDetector det;
  BarycenterInput baryinput;

  REAL8Vector *phiMean=NULL;
  REAL8 intre=0., intim=0., integral1=0.;

  REAL8 maxmismatch = 1e-100;
  REAL8 meanmismatch = 0.;

  INT4Vector *matPos=NULL; /* position of parameters in ParamData */
  UINT4 seed=0;              /* set to get seed from clock time */
  RandomParams *randomParams=NULL;
  ParamData *paramData=NULL, *randVals=NULL;

  REAL8Array *cormat=NULL, *covmat=NULL, *posdef=NULL;
  REAL8Array *chol=NULL;

  CHAR outputfile[256];

  get_input_args(&inputs, argc, argv);

  npoints = inputs.timespan/inputs.deltat;

  XLALReadTEMPOParFile( &params, inputs.parfile );
  XLALReadTEMPOParFile( &deltaparams, inputs.parfile );

  if( strstr(inputs.parfile2, ".par") != NULL )
    XLALReadTEMPOParFile( &params2, inputs.parfile2 );

  /* set up ephemerises */
  det = *XLALGetSiteInfo( inputs.det ); /* just set site as LHO */
  baryinput.site = det;
  baryinput.site.location[0] /= LAL_C_SI;
  baryinput.site.location[1] /= LAL_C_SI;
  baryinput.site.location[2] /= LAL_C_SI;

  edat = XLALMalloc(sizeof(*edat));
  (*edat).ephiles.earthEphemeris = inputs.earth;
  (*edat).ephiles.sunEphemeris = inputs.sun;
  LALInitBarycenter(&status, edat);

  /* set the position and frequency epochs if not already set */
  if(params.pepoch == 0. && params.posepoch != 0.)
    params.pepoch = params.posepoch;
  else if(params.posepoch == 0. && params.pepoch != 0.)
    params.posepoch = params.pepoch;

  if( strstr(inputs.parfile2, ".par") != NULL ){
    if(params2.pepoch == 0. && params2.posepoch != 0.)
    params2.pepoch = params2.posepoch;
  else if(params2.posepoch == 0. && params2.pepoch != 0.)
    params2.posepoch = params2.pepoch;
  }

  /* set the position and frequency epochs if not already set */
  if(deltaparams.pepoch == 0. && deltaparams.posepoch != 0.)
    deltaparams.pepoch = deltaparams.posepoch;
  else if(deltaparams.posepoch == 0. && deltaparams.pepoch != 0.)
    deltaparams.posepoch = deltaparams.pepoch;

  /* calculate the phase every minute for the mean values */
  phiMean = get_phi( inputs.start, (double)inputs.deltat, npoints, params, baryinput,
    edat);

  /* calculate phase integral over time */
  for( k=0; k<npoints-1; k++ ){
    intre += 2.;
    intim += 0.;
  }

  integral1 = intre*intre + intim*intim;

  /* see which parameters have an associated error */
  if( params.f0Err != 0. ) count++;
  if( params.f1Err != 0. ) count++;
  if( params.f2Err != 0. ) count++;
  if( params.raErr != 0. ) count++;
  if( params.decErr != 0. ) count++;
  if( params.pmraErr != 0. ) count++;
  if( params.pmdecErr != 0. ) count++;

  /* see if binary model */
  if( params.model != NULL ){
    if( params.xErr != 0. ) countBin++;
    if( params.w0Err != 0. ) countBin++;
    if( params.eErr != 0. ) countBin++;
    if( params.PbErr != 0. ) countBin++;
    if( params.T0Err != 0. ) countBin++;
    if( params.TascErr != 0. ) countBin++;
    if( params.eps1Err != 0. ) countBin++;
    if( params.eps2Err != 0. ) countBin++;
    if( params.gammaErr != 0. ) countBin++;
    if( params.wdotErr != 0. ) countBin++;
    if( params.PbdotErr != 0. ) countBin++;
    if( params.xdotErr != 0. ) countBin++;
    if( params.edotErr != 0. ) countBin++;
    if( params.sErr != 0. ) countBin++;
  }

  /* set up random parameters */
  if ( strstr(inputs.parfile2, ".par") == NULL ){
    randomParams = XLALCreateRandomParams( seed );

    paramData = XLALMalloc( MAXPARAMS*sizeof(ParamData) );

    /* get correlation matrix */
    cormat = ReadCorrelationMatrix( inputs.covfile, params, paramData );

    if( (posdef = XLALCheckPositiveDefinite( cormat )) == NULL ){
      /* set covariance matrix */
      covmat = CreateCovarianceMatrix( paramData, cormat );

      chol = CholeskyDecomp(covmat, "lower");
    }
    else{
      covmat = CreateCovarianceMatrix( paramData, posdef );

      chol = CholeskyDecomp(covmat, "lower");

      XLALDestroyREAL8Array(posdef);
    }

    /* get the posistions of the parameters within paramData */
    i=0;
    j=0;
    matPos = XLALCreateINT4Vector( MAXPARAMS );
    while(i < MAXPARAMS){
      if( paramData[i].matPos != 0 ){
        matPos->data[j] = i;
        j++;
      }
      i++;
    }
    matPos = XLALResizeINT4Vector(matPos, j);
  }

  N = inputs.iterations;

  if( strstr(inputs.parfile2, ".par") != NULL ) N = 1;

  /* calculate the max mismatch between the mean value and the mean plus the
     standard deviation */
  for( i=0; i < N; i++ ){
    REAL8Vector *phiOffset=NULL;
    REAL8 integral2=0.;

    /* generate new pulsars parameters from the pulsar parameter covariance */
    if( strstr(inputs.parfile2, ".par") == NULL ){
      randVals = MultivariateNormalDeviates( chol, paramData, randomParams );
      SetMCMCPulsarParams( &deltaparams, randVals );
    }

    /* 
    if( i==0 ){
      if(count > 0){
        if( params.f0Err != 0. ) fprintf(stderr, "Search over f0\n");
        if( params.f1Err != 0. ) fprintf(stderr, "Search over f1\n");
        if( params.f2Err != 0. ) fprintf(stderr, "Search over f2\n");
        if( params.raErr != 0. ) fprintf(stderr, "Search over ra\n");
        if( params.decErr != 0. ) fprintf(stderr, "Search over dec\n");
        if( params.pmraErr != 0. ) fprintf(stderr, "Search over pmra\n");
        if( params.pmdecErr != 0. ) fprintf(stderr, "Search over pmdec\n");
      }

      if( countBin > 0 ){
        if( params.xErr != 0. ) fprintf(stderr, "Search over x\n");
        if( params.w0Err != 0. ) fprintf(stderr, "Search over w0\n");
        if( params.eErr != 0. ) fprintf(stderr, "Search over e\n");
        if( params.PbErr != 0. ) fprintf(stderr, "Search over Pb\n");
        if( params.T0Err != 0. ) fprintf(stderr, "Search over T0\n");
        if( params.TascErr != 0. ) fprintf(stderr, "Search over Tasc\n");
        if( params.eps1Err != 0. ) fprintf(stderr, "Search over eps1\n");
        if( params.eps2Err != 0. ) fprintf(stderr, "Search over eps2\n");
        if( params.gammaErr != 0. ) fprintf(stderr, "Search over gamma\n");
        if( params.wdotErr != 0. ) fprintf(stderr, "Search over wdot\n");
        if( params.PbdotErr != 0. ) fprintf(stderr, "Search over Pbdot\n");
        if( params.xdotErr != 0. ) fprintf(stderr, "Search over xdot\n");
        if( params.edotErr != 0. ) fprintf(stderr, "Search over edot\n");
        if( params.sErr != 0. ) fprintf(stderr, "Search over sini\n");
      }
    }
    */

    /* get new phase */
    if( strstr(inputs.parfile2, ".par") == NULL ){
      phiOffset = get_phi( inputs.start, (double)inputs.deltat, npoints,
        deltaparams, baryinput, edat );
    }
    else{
      phiOffset = get_phi( inputs.start, (double)inputs.deltat, npoints,
        params2, baryinput, edat );
    }

    /* calculate the mismatch 1 - (P(params + delta))/P(params) */
    intre = 0.;
    intim = 0.;

    for( k=0; k<npoints-1; k++ ){
      REAL8 phi1=0., phi2=0.;
      phi1 = LAL_TWOPI*fmod(phiOffset->data[k]-phiMean->data[k], 1.);
      phi2 = LAL_TWOPI*fmod(phiOffset->data[k+1]-phiMean->data[k+1], 1.);
      intre += (cos(phi2) + cos(phi1));
      intim += (sin(phi2) + sin(phi1));
    }

    integral2 = intre*intre + intim*intim;  /* square value for power */

    /* get the mean mismatch */
    meanmismatch += fabs(1. - integral2/integral1);

    /* work out mismatch */
    if( fabs(1. - integral2/integral1) > fabs(maxmismatch) )
      maxmismatch = 1. - integral2/integral1;

    XLALDestroyREAL8Vector( phiOffset );
  }

  meanmismatch /= (double)N;

  if( strstr(inputs.parfile2, ".par") == NULL ){
    sprintf(outputfile, "%s/mismatch_%s", inputs.outputdir, inputs.pulsar);
    fp = fopen(outputfile, "w");

    fprintf(stderr, "Maximum mismatch for %d templates drawn randomly from the\
 given parameter and pulsar %s\n%le\n", N, inputs.pulsar, maxmismatch);

    fprintf(stderr, "Mean mismatch is %le\n", meanmismatch);

    fprintf(fp, "%% Maximum and mean mismatch for %d templates drawn randomly from the \
given parameter and pulsar %s\n%lf\t%lf\n", N, inputs.pulsar, maxmismatch, meanmismatch);

    fclose(fp);
  }
  else{
    fprintf(stderr, "The mismatch between the two par file for pulsar %s is %le\n", 
      inputs.pulsar, maxmismatch);
  }

  XLALFree( paramData );
  XLALFree( randVals );
  XLALDestroyREAL8Array( chol );

  XLALDestroyINT4Vector( matPos );

  XLALDestroyREAL8Vector( phiMean );

  return 0;
}

void get_input_args(InputParams *inputParams, INT4 argc, CHAR *argv[]){
  struct option long_options[] =
  {
    { "help",           no_argument,       0, 'h' },
    { "detector",      required_argument, 0, 'D' },
    { "pulsar",         required_argument, 0, 'p' },
    { "par-file",       required_argument, 0, 'P' },
    { "par-file2",      required_argument, 0, 'b' },
    { "covariance",      required_argument, 0, 'c' },
    { "start",     required_argument, 0, 's' },
    { "timespan",          required_argument, 0, 't' },
    { "deltat", required_argument, 0, 'd' },
    { "iterations", required_argument, 0, 'i' },
    { "output-dir", required_argument, 0, 'o' },
    { "earth", required_argument, 0, 'e' },
    { "sun", required_argument, 0, 'S' },
    { 0, 0, 0, 0}
  };

  CHAR args[] = "hd:p:P:b:c:s:t:d:i:o:e:S:" ;
  CHAR *program = argv[0];

  while( 1 ){
    INT4 option_index = 0;
    INT4 c;

    c = getopt_long_only( argc, argv, args, long_options, &option_index );
    if( c == -1 ) /* end of options */
      break;

    switch( c ){
      case 0:
        if( long_options[option_index].flag )
          break;
        else
          fprintf(stderr, "Error passing option %s with argument %s\n",
            long_options[option_index].name, optarg);
      case 'h': /* help message */
        fprintf(stderr, USAGE, program);
        exit(0);
      case 'D': 
        sprintf(inputParams->det, "%s", optarg);
        break;
      case 'p':
        sprintf(inputParams->pulsar, "%s", optarg);
        break;
      case 'P':
        sprintf(inputParams->parfile, "%s", optarg);
        break;
      case 'b':
        sprintf(inputParams->parfile2, "%s", optarg);
        break;
      case 'c':
        sprintf(inputParams->covfile, "%s", optarg);
        break;
      case 's':
        inputParams->start = atoi(optarg);
        break;
      case 't':
        inputParams->timespan = atoi(optarg);
        break;
      case 'd':
        inputParams->deltat = atoi(optarg);
        break;
      case 'i':
        inputParams->iterations = atoi(optarg);
        break;
      case 'o':
        sprintf(inputParams->outputdir, "%s", optarg);
        break;
      case 'e':
        sprintf(inputParams->earth, "%s", optarg);
        break;
      case 'S':
        sprintf(inputParams->sun, "%s", optarg);
        break;
      case '?':
        fprintf(stderr, "Unknown error while parsing options\n");
      default:
        fprintf(stderr, "Unknown error while parsing options\n");
    }
  }
}

/* function to return a vector of the pulsar phase for each data point */
REAL8Vector *get_phi( double start, double deltaT, int npoints,
 BinaryPulsarParams params, BarycenterInput bary, EphemerisData *edat ){
  static LALStatus status;

  INT4 i=0;

  REAL8 T0=0., DT=0., DTplus=0., deltat=0., deltat2=0.;
  REAL8 interptime = 1800.; /* calulate every 30 mins (1800 secs) */
  REAL8 time=0.;

  EarthState earth, earth2;
  EmissionTime emit, emit2;
  REAL8 emitdt=0.;

  BinaryPulsarInput binput;
  BinaryPulsarOutput boutput;

  REAL8Vector *phis=NULL;

  /* if edat is NULL then return a NULL poniter */
  if( edat == NULL)
    return NULL;

  /* allocate memory for phases */
  phis = XLALCreateREAL8Vector( npoints );

  bary.dInv = 0.;

  for( i=0; i<npoints; i++){
    T0 = params.pepoch;
    time = start + deltaT*(double)i;

    DT = time - T0;

    if(time <= 820108813)
      (*edat).leap = 13;
    else if(time <= 914803214)
      (*edat).leap = 14;
    else
      (*edat).leap = 15;

    /* only call the barycentring routines every 30 minutes, otherwise just
       linearly interpolate between them */
    if( i==0 || DT > DTplus ){
      bary.tgps.gpsSeconds = (UINT8)floor(time);
      bary.tgps.gpsNanoSeconds = (UINT8)floor((fmod(time,1.)*1e9));

      bary.delta = params.dec + DT*params.pmdec;
      bary.alpha = params.ra + DT*params.pmra/cos(bary.delta);

      /* call barycentring routines */
      LALBarycenterEarth(&status, &earth, &bary.tgps, edat);
      LALBarycenter(&status, &emit, &bary, &earth);

      /* add 30 minutes (1800secs) to the time */
      DTplus = DT + interptime;
      bary.tgps.gpsSeconds = (UINT8)floor(time+interptime);
      bary.tgps.gpsNanoSeconds = (UINT8)floor((fmod(time+interptime,1.)*1e9));

      /* No point in updating the positions as difference will be tiny */
      LALBarycenterEarth(&status, &earth2, &bary.tgps, edat);
      LALBarycenter(&status, &emit2, &bary, &earth2);
    }

    /* linearly interpolate to get emitdt */
    emitdt = emit.deltaT + (DT - (DTplus - interptime)) *
      (emit2.deltaT - emit.deltaT)/interptime;

    /* check if need to perform binary barycentring */
    if( params.model != NULL ){
      binput.tb = time + emitdt;

      XLALBinaryPulsarDeltaT( &boutput, &binput, &params );

      deltat = DT + emitdt + boutput.deltaT;
    }
    else
      deltat = DT + emitdt;

    /* work out phase */
    deltat2 = deltat*deltat;
    phis->data[i] = 2.*(params.f0*deltat + 0.5*params.f1*deltat2
      + SIXTH*params.f2*deltat*deltat2
      + TWENTYFOURTH*params.f3*deltat2*deltat2
      + (1./120.)*params.f4*deltat2*deltat2 
      + (1./720.)*params.f5*deltat2*deltat2*deltat);
  }

  return phis;
}

void SetMCMCPulsarParams( BinaryPulsarParams *pulsarParams, ParamData *data ){
  /* go through params and set the next step in the MCMC */
  pulsarParams->f0 = data[0].val;
  pulsarParams->f1 = data[1].val;
  pulsarParams->f2 = data[2].val;
  pulsarParams->dec = data[3].val;
  pulsarParams->ra = data[4].val;
  pulsarParams->pmdec = data[5].val;
  pulsarParams->pmra = data[6].val;
  pulsarParams->x = data[7].val;
  pulsarParams->e = data[8].val;
  pulsarParams->T0 = data[9].val;
  pulsarParams->Pb = data[10].val;
  pulsarParams->w0 = data[11].val;
  pulsarParams->wdot = data[12].val;
  pulsarParams->gamma = data[13].val;
  pulsarParams->Pbdot = data[14].val;
  pulsarParams->s = data[15].val;
  pulsarParams->M = data[16].val;
  pulsarParams->m2 = data[17].val;
  pulsarParams->dth = data[18].val;
  pulsarParams->xdot = data[19].val;
  pulsarParams->edot = data[20].val;
  pulsarParams->x2 = data[21].val;
  pulsarParams->e2 = data[22].val;
  pulsarParams->T02 = data[23].val;
  pulsarParams->Pb2 = data[24].val;
  pulsarParams->w02 = data[25].val;
  pulsarParams->x3 = data[26].val;
  pulsarParams->e3 = data[27].val;
  pulsarParams->T03 = data[28].val;
  pulsarParams->Pb3 = data[29].val;
  pulsarParams->w03 = data[30].val;
  pulsarParams->xpbdot = data[31].val;
  pulsarParams->f3 = data[32].val;
  pulsarParams->f4 = data[33].val;
  pulsarParams->f5 = data[34].val;

  if( pulsarParams->model != NULL && !strcmp(pulsarParams->model, "ELL1") ){
    pulsarParams->eps1 = data[8].val;
    pulsarParams->eps2 = data[11].val;
    pulsarParams->e = 0.;
    pulsarParams->Tasc = data[9].val;
    pulsarParams->eps1dot = data[12].val;
    pulsarParams->eps2dot = data[20].val;
  }
}

/* this function perform Cholesky decomposition on M and outputs the
   lower, or upper triagular matrix depending if uOrl is set to "upper" or
   "lower" - if nothing is specified then the default is lower
   This is pretty much copied from the GSL function gsl_linalg_cholesky_decomp
   although this works with floats rather than doubles
*/
REAL8Array *CholeskyDecomp(REAL8Array *M, CHAR* uOrl){
  INT4 i=0, j=0, k=0;

  REAL8 A_00=0., L_00=0.;

  REAL8Array *A=NULL;

  INT4 length=0;

  /* fprintf(stderr, "Performing Cholesky decomposition of matrix\n"); */

  /* check dimensions are equal */
  if( M->dimLength->data[0] != M->dimLength->data[1] ){
    fprintf(stderr, "Error... input matrix has unequal dimensions!\n");
    exit(0);
  }

  length = M->dimLength->data[0];

  /* allocate memory */
  A = XLALCreateREAL8Array( M->dimLength );

  if(M == NULL || A == NULL){
    fprintf(stderr, "Error... input or output matrix is NULL!\n");
    exit(0);
  }

  /* initialise L be same as input matrix M */
  for(i=0; i < length; i++)
    for(j=0; j < length; j++)
      A->data[i*length + j] = M->data[i*length + j];

  A_00 = A->data[0];
  L_00 = sqrt(A_00);

  if( A_00 <= 0 )
    fprintf(stderr, "Error... matrix must be positive definite!\n");

  A->data[0] = L_00;

  if( length > 1 ){
    REAL8 A_10 = A->data[1*length + 0];
    REAL8 A_11 = A->data[1*length + 1];

    REAL8 L_10 = A_10/L_00;
    REAL8 diag = A_11 - L_10*L_10;
    REAL8 L_11 = sqrt(diag);

    if( diag <= 0 ){
      fprintf(stderr, "Error... input matrix is not pos def!\n");
      exit(0);
    }

    A->data[1*length + 0] = L_10;
    A->data[1*length + 1] = L_11;
  }

  for( k=2; k<length; k++ ){
    REAL8 A_kk = A->data[k*length + k];

    for( i=0; i<k; i++ ){
      REAL8 sum = 0.;

      REAL8 A_ki = A->data[k*length + i];
      REAL8 A_ii = A->data[i*length + i];

      REAL8 ci[length];
      REAL8 ck[length];

      for( j=0; j<length; j++ ){
        ci[j] = A->data[i*length + j];
        ck[j] = A->data[k*length + j];
      }

      if( i>0 ){
        for( j=0; j<i; j++ )
          sum += ci[j]*ck[j];
      }

      A_ki = (A_ki - sum) / A_ii;
      A->data[k*length + i] = A_ki;
    }

    {
      REAL8 sum = 0.;
      REAL8 diag = 0.;

      for( j=0; j<k; j++ ){
        sum += A->data[k*length + j] * A->data[k*length + j];
      }

      diag = A_kk - sum;

      /* check if the diag value is negative, but also not close to the minimum
         resolvable difference between to REAL8 numbers - if it is negative
         and also close to this value set it to +LAL_REAL8_EPS (see
         LALConstants.h), as it could be anywhere inbetween LAL_REAL8_MIN and
         LAL_REAL8_EPS */
      /* these problems are caused by the fact the when computing the eigen
         values/vectors to determine if the matrix is positive definite the
         process uses iterative methods to check on the convergence of values
         and these will only be good down to the precision of REAL8s */
      if( diag <= 0. && fabs(diag) >= LAL_REAL8_EPS && k != length-1 ){
        fprintf(stderr, "Error... input matrix is not pos def!\n");
        exit(0);
      }
      else if( diag <= 0. && fabs(diag) <= LAL_REAL8_EPS ){
        diag = fabs(diag);
        /* diag = LAL_REAL8_EPS; */
        /* diag = 0.; */
      }
      else if( diag <= 0. && fabs(diag) >= LAL_REAL8_EPS && k == length-1 ){
        /* this is a kludge as a lot of the matricies seem to have entries
           there m(length, length) diagonal value as small but less than zero,
           so I'll just set it to zero manually */
        diag = 0.;
      }

      A->data[k*length + k] = sqrt(diag);

    }
  }

  /* set upper triangular matrix to zeros - for lower value */
  for(i=0; i<length; i++)
    for(j=i+1; j<length; j++)
      A->data[i*length + j] = 0.;

  /* check if the upper triangle is wanted - if so perform transpose */
  if(strstr(uOrl, "upper")!=NULL){
    REAL8 tempdata = 0.;

    /* perform transpose */
    for(j=0; j<length-1; j++){
      for(i=j; i<length; i++){
        tempdata = A->data[j*length + i];
        A->data[j*length + i] = A->data[i*length + j];
        A->data[i*length + j] = tempdata;
      }
    }
  }

  /* fprintf(stderr, "\nCholesky decomposed matrix:\n");
  for(i=0; i<length; i++){
     for(j=0; j<length; j++)
       fprintf(stderr, "%.2e  ", A->data[i*length + j]);
     fprintf(stderr, "\n");
  } */

  return A;
}

/* this function will draw a set of random numbers from a multivariate Gaussian
distribution, with a cholesky decomposed covariance matrix given by cholmat and
parameter mean values */
ParamData *MultivariateNormalDeviates( REAL8Array *cholmat, ParamData *data,
  RandomParams *randomParams ){
  REAL4Vector *randNum=NULL;

  ParamData *deviates=NULL;

  INT4 dim=cholmat->dimLength->data[0]; /* covariance matrix dimensions */

  INT4 i=0, j=0;
  REAL8Vector *Z=NULL;

  /* check dimensions of covariance matrix and mean vector length are equal */
  if( cholmat->dimLength->data[0] != cholmat->dimLength->data[1] ){
    fprintf(stderr, "Error... wrong number of dimensions in input matrix!\n");
    exit(0);
  }

  deviates = XLALMalloc(MAXPARAMS*sizeof(ParamData));

  /* create a vector of random Gaussian numbers */
  randNum = XLALCreateREAL4Vector( dim );

  XLALNormalDeviates( randNum, randomParams );

  /* multiply L by randNum */
  Z = XLALCreateREAL8Vector( dim );
  for(i=0;i<dim;i++)
    Z->data[i] = 0.;

  for(i=0;i<dim;i++)
    for(j=0;j<dim;j++)
      Z->data[i] += cholmat->data[i*dim +j]*randNum->data[j];

  /* get the output random deviates by doing the mean plus Z */
  j=0;
  for(i=0;i<MAXPARAMS;i++){
    deviates[i].name = data[i].name;
    deviates[i].sigma = data[i].sigma;
    deviates[i].matPos = data[i].matPos;
    if( data[i].matPos != 0 )
      deviates[i].val = data[i].val + Z->data[data[i].matPos-1];
    else
      deviates[i].val = data[i].val;
  }

  XLALDestroyREAL4Vector( randNum );
  XLALDestroyREAL8Vector( Z );

  return deviates;
}

/* I need to define a standard set of positions in which various pulsar
   parameters will sit within the internal correlation matrix - this will as
   far as possible following the standard in the matrix files I have from
   Michael Kramer
*/
/* read in the correlation matrix */
REAL8Array *ReadCorrelationMatrix( CHAR *matrixFile, 
  BinaryPulsarParams params, ParamData *data ){
  FILE *fp=NULL;

  CHAR matrixParams[MAXPARAMS][6]; /* parameters in the correlation matrix */
  CHAR paramTmp[256];

  INT4 numParams=0, i=0, j=0, k=0, n=0;

  CHAR tmpStr[256], tmpStr2[256];

  INT4 arraySize=0;

  REAL8Array *corMat=NULL;
  UINT4Vector *matdims=NULL;

  INT4 DMpos=0, DM1pos=0; /* position of dispersion measure in matrix */
  INT4 numDM=0;
  REAL8 corTemp=0., junk=0.;

  ParamData paramData[]=
  {
    { "f0",  0., 0., 0 },{ "f1",  0., 0., 0 },{ "f2",  0., 0., 0 },
    { "Dec", 0., 0., 0 },{ "RA",  0., 0., 0 },{ "pmdc",0., 0., 0 },
    { "pmra",0., 0., 0 },{ "x",   0., 0., 0 },{ "e",   0., 0., 0 },
    { "T0",  0., 0., 0 },{ "Pb",  0., 0., 0 },{ "Om",  0., 0., 0 },
    { "Omdt",0., 0., 0 },{ "gamma",0., 0., 0 },{ "Pbdt",0.,0., 0 },
    { "s",   0., 0., 0 },{ "M",   0., 0., 0 },{ "m2",  0., 0., 0 },
    { "dth", 0., 0., 0 },{ "xdot",0., 0., 0 },{ "edot",0., 0., 0 },
    { "x2",  0., 0., 0 },{ "e2",  0., 0., 0 },{ "T02", 0., 0., 0 },
    { "Pb2", 0., 0., 0 },{ "Om2", 0., 0., 0 },{ "x3",  0., 0., 0 },
    { "e3",  0., 0., 0 },{ "T03", 0., 0., 0 },{ "Pb3", 0., 0., 0 },
    { "Om3", 0., 0., 0 },{ "Xpbd",0., 0., 0 },{ "f3",  0., 0., 0 },
    { "f4",  0., 0., 0 },{ "f5",  0., 0., 0 }
  };

  /* set the values - put the parameter errors at twice those given in the .par
files */
  paramData[0].val = params.f0;      paramData[0].sigma = params.f0Err;
  paramData[1].val = params.f1;      paramData[1].sigma = params.f1Err;
  paramData[2].val = params.f2;      paramData[2].sigma = params.f2Err;
  paramData[3].val = params.dec;     paramData[3].sigma = params.decErr;
  paramData[4].val = params.ra;      paramData[4].sigma = params.raErr;
  paramData[5].val = params.pmdec;   paramData[5].sigma = params.pmdecErr;
  paramData[6].val = params.pmra;    paramData[6].sigma = params.pmraErr;
  paramData[7].val = params.x;       paramData[7].sigma = params.xErr;
  paramData[8].val = params.e;       paramData[8].sigma = params.eErr;
  paramData[9].val = params.T0;      paramData[9].sigma = params.T0Err;
  paramData[10].val = params.Pb;     paramData[10].sigma = params.PbErr;
  paramData[11].val = params.w0;     paramData[11].sigma = params.w0Err;
  paramData[12].val = params.wdot;   paramData[12].sigma = params.wdotErr;
  paramData[13].val = params.gamma;  paramData[13].sigma = params.gammaErr;
  paramData[14].val = params.Pbdot;  paramData[14].sigma = params.PbdotErr;
  paramData[15].val = params.s;      paramData[15].sigma = params.sErr;
  paramData[16].val = params.M;      paramData[16].sigma = params.MErr;
  paramData[17].val = params.m2;     paramData[17].sigma = params.m2Err;
  paramData[18].val = params.dth;    paramData[18].sigma = params.dthErr;
  paramData[19].val = params.xdot;   paramData[19].sigma = params.xdotErr;
  paramData[20].val = params.edot;   paramData[20].sigma = params.edotErr;
  paramData[21].val = params.x2;     paramData[21].sigma = params.x2Err;
  paramData[22].val = params.e2;     paramData[22].sigma = params.e2Err;
  paramData[23].val = params.T02;    paramData[23].sigma = params.T02Err;
  paramData[24].val = params.Pb2;    paramData[24].sigma = params.Pb2Err;
  paramData[25].val = params.w02;    paramData[25].sigma = params.w02Err;
  paramData[26].val = params.x3;     paramData[26].sigma = params.x3Err;
  paramData[27].val = params.e3;     paramData[27].sigma = params.e3Err;
  paramData[28].val = params.T03;    paramData[28].sigma = params.T03Err;
  paramData[29].val = params.Pb3;    paramData[29].sigma = params.Pb3Err;
  paramData[30].val = params.w03;    paramData[30].sigma = params.w03Err;
  paramData[31].val = params.xpbdot; paramData[31].sigma = params.xpbdotErr;
  paramData[32].val = params.f3;     paramData[32].sigma = params.f3Err;
  paramData[33].val = params.f4;     paramData[33].sigma = params.f4Err;
  paramData[34].val = params.f5;     paramData[34].sigma = params.f5Err;

  arraySize = MAXPARAMS;

  if( params.model != NULL && !strcmp(params.model, "ELL1") ){
    paramData[8].name = "eps1";
    paramData[8].val = params.eps1;
    paramData[8].sigma = params.eps1Err;

    paramData[11].name = "eps2";
    paramData[11].val = params.eps2;
    paramData[11].sigma = params.eps2Err;

    paramData[9].name = "Tasc";
    paramData[9].val = params.Tasc;
    paramData[9].sigma = params.TascErr;

    paramData[12].name = "e1dt";
    paramData[12].val = params.eps1dot;
    paramData[12].sigma = params.eps1dotErr;

    paramData[20].name = "e2dt";
    paramData[20].val = params.eps2dot;
    paramData[20].sigma = params.eps2dotErr;
  }

  /* read in data from correlation matrix file */
  if((fp = fopen(matrixFile, "r")) == NULL){
    fprintf(stderr, "No correlation matrix file" );
    return NULL;
  }

  /* read in the first line of the matrix file */
  while(fscanf(fp, "%s", paramTmp)){
    if(strchr(paramTmp, '-') != NULL)
      break;

    if(feof(fp)){
      fprintf(stderr, "Error... I've reached the end of the file without \
reading any correlation data!");
      fclose(fp);
      exit(0);
    }

    sprintf(matrixParams[numParams], "%s", paramTmp);
    numParams++;

    /* check if parameter is actually for a dispersion measure (ignore if so) */
    if(!strcmp(paramTmp, "DM")){
      numParams--;
      DMpos = i;
      numDM++;
    }
    if(!strcmp(paramTmp, "DM1")){
      numParams--;
      DM1pos = i;
      numDM++;
    }

    i++;
  };

  if(numParams > arraySize){
    fprintf(stderr, "More parameters in matrix file than there should be!\n");
    exit(0);
  }

  matdims = XLALCreateUINT4Vector( 2 );
  matdims->data[0] = numParams;
  matdims->data[1] = numParams;

  corMat = XLALCreateREAL8Array( matdims );

  /* find positions of each parameter */
  /* the strings that represent parameters in a matrix are given in the paraem
     variable in the tempo code mxprt.f */
  /* read in matrix */
  k=0;
  for(i=0;i<numParams+numDM;i++){
    n=0;
    fscanf(fp, "%s%s", tmpStr, tmpStr2);

    /* if its a dispersion measure then just skip the line */
    if( (DMpos != 0 && i == DMpos) || (DM1pos != 0 && i == DM1pos) ){
      fscanf(fp, "%*[^\n]");
      k--;
      continue;
    }

    for(j=0;j<i+1;j++){
      if( (DMpos != 0 && j == DMpos) || (DM1pos != 0 && j == DM1pos) ){
        fscanf(fp, "%lf", &junk);
        n--;
        continue;
      }

      fscanf(fp, "%lf", &corTemp);
      
      /* if covariance equals 1 set as 0.9999999, because values of 1
         can cause problems of giving singular matrices */
      if( (n != k) && (corTemp == 1.) )
        corTemp = 0.9999999;
      else if( (n != k) && (corTemp == -1.) )
        corTemp = -0.9999999;

      corMat->data[k*corMat->dimLength->data[0] + n] = corTemp;
      if(n != k)
        corMat->data[n*corMat->dimLength->data[0] + k] = corTemp;
      n++;
    }

    /* send an error if we hit the end of the file */
    if(feof(fp)){
      fprintf(stderr, "Error reading in matrix - hit end of file!\n");
      exit(0);
    }

    k++;
  }

  fclose(fp);

  /* give the correlation matrix positions of the parameters */
  for(i=1;i<numParams+1;i++){
    for(j=0;j<arraySize;j++){
      if(!strcmp(matrixParams[i-1], paramData[j].name)){
        paramData[j].matPos = i;

        break;
      }
    }
  }

  /* pass the parameter data to be output */
  memcpy(data, paramData, sizeof(paramData));

  XLALDestroyUINT4Vector( matdims );

  return corMat;
}

/* turn the input correlation matrix into a covariance matrix */
REAL8Array *CreateCovarianceMatrix( ParamData *data, REAL8Array *corMat ){
  REAL8Array *covMat=NULL;
  INT4 i=0, j=0;

  covMat = XLALCreateREAL8Array( corMat->dimLength );

  /* convert correlation matrix into a covariance matrix */
  for(i=0;i<MAXPARAMS;i++){
    if( data[i].matPos != 0 ){
      for(j=0;j<MAXPARAMS;j++){
        if( data[j].matPos != 0 ){
          covMat->data[(data[i].matPos-1)*covMat->dimLength->data[0] +
            data[j].matPos-1] =
          corMat->data[(data[i].matPos-1)*corMat->dimLength->data[0] +
             data[j].matPos-1] * data[i].sigma * data[j].sigma;
        }
      }
    }
  }

  return covMat;
}

REAL8Array *XLALCheckPositiveDefinite( REAL8Array *matrix ){
  static LALStatus status;

  REAL8Vector *eigenval=NULL;
  REAL8Array *eigenvec=NULL;

  REAL8Array *posdef=NULL;

  INT4 i=0, j=0;

  /* copy input array into eigenvec as this gets converted by function */
  eigenvec = XLALCreateREAL8Array( matrix->dimLength );

  for( i=0; i<(INT4)eigenvec->dimLength->data[0]; i++ ){
    for( j=0; j<(INT4)eigenvec->dimLength->data[1]; j++ ){
      eigenvec->data[i*eigenvec->dimLength->data[0] + j] =
        matrix->data[i*matrix->dimLength->data[0] + j];
    }
  }

  eigenval = XLALCreateREAL8Vector( matrix->dimLength->data[0] );

  /* calculate the eigen values and vectors */
  LALDSymmetricEigenVectors( &status, eigenval, eigenvec );

  for( i=0; i<(INT4)matrix->dimLength->data[0]; i++ ){
    /* first check if any eigen values are zero and if so convert to positive
       definite matrix */
    if( eigenval->data[i] < 0. && fabs(eigenval->data[i]) > 10.*LAL_REAL8_EPS ){
      fprintf(stderr, "Eigenvalue is negative. Non-postive definite matrix!\n");
      posdef = XLALConvertToPositiveDefinite( matrix );
      break;
    }
  }

  /* if matrix is positive definite return it i.e. posdef hasn't been set */
  if( posdef == NULL ){
    XLALDestroyREAL8Array( eigenvec );
    XLALDestroyREAL8Vector( eigenval );
    return NULL;
  }

  /* re-check new matrix for positive definiteness, but be aware of values
     close to the precision of REAL8 numbers */
  for( i=0; i<(INT4)eigenvec->dimLength->data[0]; i++ ){
    for( j=0; j<(INT4)eigenvec->dimLength->data[1]; j++ ){
      eigenvec->data[i*eigenvec->dimLength->data[0] + j] =
        posdef->data[i*posdef->dimLength->data[0] + j];
    }
    eigenval->data[i] = 0.;
  }

  LALDSymmetricEigenVectors( &status, eigenval, eigenvec );

  for( i=0; i<(INT4)matrix->dimLength->data[0]; i++ ){
    if( eigenval->data[i] < 0. && fabs(eigenval->data[i]) > 10.*LAL_REAL8_EPS){
      fprintf(stderr, "ABORT! Eigenvalue is negative. Non-postive definite \
matrix!\n");
      exit(0);
    }
  }

  XLALDestroyREAL8Array( eigenvec );
  XLALDestroyREAL8Vector( eigenval );

  return posdef;
}

/* this function takes a matrix that isn't positive definite and converts it
into a postive definite matrix using the method (number 2) of Rebonato and
Jackel (see their paper at
http://www.riccardorebonato.co.uk/papers/ValCorMat.pdf */
REAL8Array *XLALConvertToPositiveDefinite( REAL8Array *nonposdef ){
  static LALStatus status;

  REAL8Vector *eigenval=NULL;
  REAL8Array *eigenvec=NULL;

  REAL8Array *posdef = NULL; /* output positive definite matrix */

  REAL8Array *Lprime=NULL;
  REAL8Array *B=NULL, *Bprime=NULL, *Btrans=NULL;
  REAL8Array *T=NULL;

  REAL8 Tval=0.;

  INT4 i=0, j=0, length=0;

  fprintf(stderr, "Converting to positive definite matrix\n");

  /* check that matrix is square */
  if( nonposdef->dimLength->data[0] != nonposdef->dimLength->data[1] ){
    fprintf(stderr, "Error... matrix must be square!\n");
    exit(0);
  }
  else
    length = nonposdef->dimLength->data[0];

  Lprime = XLALCreateREAL8Array( nonposdef->dimLength );
  T = XLALCreateREAL8Array( nonposdef->dimLength );

  /* copy input array into eigenvec as this gets converted by function */
  eigenvec = XLALCreateREAL8Array( nonposdef->dimLength );

  for( i=0; i<length; i++ ){
    for( j=0; j<length; j++ ){
      eigenvec->data[i*length + j] = nonposdef->data[i*length + j];

      /* initialise Lprime and T to zeros */
      Lprime->data[i*length + j] = 0.;
      T->data[i*length + j] = 0.;
    }
  }

  eigenval = XLALCreateREAL8Vector( length );

  /* calculate the eigen values and vectors */
  LALDSymmetricEigenVectors( &status, eigenval, eigenvec );

  /* if eigen value is > 0 set Lprime to that value i.e. have eigen values of 
     zero if eigen value is negative */
  for( i=0; i<length; i++ )
    if( eigenval->data[i] > 0. )
      Lprime->data[i*length + i] = eigenval->data[i];

  /* compute scaling matrix T */
  for( i=0; i<length; i++ ){
    Tval = 0.;
    for( j=0; j<length; j++ ){
      Tval += eigenvec->data[i*length + j] * eigenvec->data[i*length + j] *
        Lprime->data[j*length + j];
    }

    Tval = 1./Tval;

    /* really we just want the sqrt of T */
    T->data[i*length + i] = sqrt(Tval);
  }

  /* convert Lprime to sqrt(lambdaprime) */
  for( i=0; i<length; i++ ){
    REAL8 tempL = Lprime->data[i*length + i];

    Lprime->data[i*length + i] = sqrt(tempL);
  }

  /* Bprime = S*sqrt(lambdaprime); */
  Bprime = XLALCreateREAL8Array( nonposdef->dimLength );
  LALDMatrixMultiply( &status, Bprime, eigenvec, Lprime );

  /* B = sqrt(T)*Bprime */
  B = XLALCreateREAL8Array( nonposdef->dimLength );
  LALDMatrixMultiply( &status, B, T, Bprime );

  /* transpose(B) */
  Btrans = XLALCreateREAL8Array( nonposdef->dimLength );
  LALDMatrixTranspose( &status, Btrans, B );

  /* posdef = B*transpose(B) this is our new positive definite matrix */
  posdef = XLALCreateREAL8Array( nonposdef->dimLength );
  LALDMatrixMultiply( &status, posdef, B, Btrans );

  XLALDestroyREAL8Array( eigenvec );
  XLALDestroyREAL8Vector( eigenval );
  XLALDestroyREAL8Array( Lprime );
  XLALDestroyREAL8Array( T );
  XLALDestroyREAL8Array( Bprime );
  XLALDestroyREAL8Array( B );
  XLALDestroyREAL8Array( Btrans );

  /* check that the difference between new and old values are greater than the
     maximum precision between REAL8 values (LAL_REAL8_EPS) - if not use
     original value */
  for( i=0; i<length; i++ ){
    for( j=0; j<length; j++ ){
      if( fabs(posdef->data[i*length + j] -
            nonposdef->data[i*length + j]) <= LAL_REAL8_EPS ){
        posdef->data[i*length + j] = nonposdef->data[i*length + j];
      }
    }
  }

  return posdef;
}
