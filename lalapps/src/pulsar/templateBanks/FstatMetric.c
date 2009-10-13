/*
 * Copyright (C) 2006 Reinhard Prix
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
 * \file 
 *
 * \author{Reinhard Prix}
 *
 * Standalone code to calculated the Fstat-metrics and mismatches
 *
 * Revision: $Id$
 *           
 */

/* ---------- includes ---------- */
#include <math.h>


#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


#include <lal/UserInput.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/LALInitBarycenter.h>
#include <lal/AVFactories.h>
#include <lal/SkyCoordinates.h>
#include <lal/ComputeFstat.h>
#include <lal/PulsarTimes.h>
#include <lal/SFTutils.h>

#include <lal/FlatPulsarMetric.h>

#include <lalapps.h>
#include <lal/ComputeFstat.h>
#include <lal/LogPrintf.h>


#include <lal/lalGitID.h>
#include <lalappsGitID.h>

RCSID ("$Id$");

/* ---------- Error codes and messages ---------- */
#define FSTATMETRIC_EMEM 	1
#define FSTATMETRIC_EINPUT	2
#define FSTATMETRIC_ENULL	3
#define FSTATMETRIC_ENONULL	4
#define FSTATMETRIC_EFILE	5


#define FSTATMETRIC_MSGEMEM	"Out of memory."
#define FSTATMETRIC_MSGEINPUT	"Invalid input."
#define FSTATMETRIC_MSGENULL	"Illegal NULL input"
#define FSTATMETRIC_MSGENONULL	"Output field not initialized to NULL"
#define FSTATMETRIC_MSGEFILE	"File IO error"

/* ----- compile switches ----- */
/* uncomment the following to turn off range-checking in GSL vector-functions */
/* #define GSL_RANGE_CHECK_OFF 1*/

#define METRIC_FORMAT	"%.15g"		/* fprintf-format used for printing metric components */

/*---------- local defines ---------- */
#define TRUE 		(1==1)
#define FALSE 		(1==0)

#define OneBillion 	1.0e9

#define NUM_SPINS 	2
#define METRIC_DIM 	2 + NUM_SPINS

#define ORB_V0		1e-4

/* ----- Macros ----- */
/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

/** Simple Euklidean scalar product for two 3-dim vectors in cartesian coords */
#define SCALAR(u,v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])

/** copy 3 components of Euklidean vector */
#define COPY_VECT(dst,src) do { (dst)[0] = (src)[0]; (dst)[1] = (src)[1]; (dst)[2] = (src)[2]; } while(0)

#define SQ(x) ((x) * (x))

/* ---------- local types ---------- */
typedef enum {
  UNITS_SI,		/**< SI-units: use seconds for 't' and for 'rX/c' */
  UNITS_NATURAL,	/**< "natural units": use Tobs for 't' and 'AU/c' for 'rX/c' */
  UNITS_LAST
} UnitsType_t;

typedef enum {
  PHASE_NONE = -1,
  PHASE_FULL = 0,
  PHASE_ORBITAL,
  PHASE_SPIN,
  PHASE_PTOLE,
  PHASE_LAST
} PhaseType_t;

typedef enum {
  METRIC_ALL = 0,
  METRIC_FSTAT,
  METRIC_PHASE,
  METRIC_ORBITAL,
  METRIC_PTOLE,
  METRIC_FLAT,
  METRIC_LAST
} MetricType_t;

/** a 'point' in the "Doppler parameter space" {alpha, delta, fkdot } */
typedef struct
{
  SkyPosition skypos;
  REAL8Vector *fkdot;
} DopplerPoint;

typedef struct
{
  REAL8Vector *dAlpha;
  REAL8Vector *dDelta;
  REAL8Vector *dFreq;
  REAL8Vector *df1dot;
} PhaseDerivs;

typedef struct
{
  UINT4 length;		/**< number of IFOs */
  PhaseDerivs **data;	/**< phase-derivs array */
} MultiPhaseDerivs;

typedef struct {
  REAL8 pos[3];
  REAL8 vel[3];
} PosVel_t;

typedef struct
{
  EphemerisData *edat;		/**< ephemeris data (from LALInitBarycenter()) */
  LIGOTimeGPS startTime;	/**< start time of observation */
  LIGOTimeGPS refTime;		/**< reference time for spin-parameters  */
  DopplerPoint dopplerPoint;	/**< sky-position and spins */
  DopplerPoint offsetUnits;	/**< (natural)units for skypos and spins */
  gsl_vector *dopplerOffset;	/**< offset-vector from signal-location */
  MultiDetectorStateSeries *multiDetStates;/**< pos, vel and LMSTs for detector at times t_i */
  MultiAMCoeffs *multiAMcoe;         	/**< Amplitude Modulation coefficients a,b(t)*/
  MultiPhaseDerivs *multidPhi;		/**< Phase-derivatives d_i phi(t) */
  REAL8Vector *GLweights;		/**< Gauss-Legendre Integration-weights */
  MultiNoiseWeights *multiNoiseWeights;	/**< noise-weights for the different IFOs */
  REAL8 Al1, Al2, Al3;			/**< amplitude factors alpha1, alpha2, alpha3 */
  REAL8 Ad, Bd, Cd, Dd;
} ConfigVariables;


typedef struct
{
  BOOLEAN help;

  LALStringVector* IFOs;	/**< list of detector-names "H1,H2,L1,.." or single detector*/
  LALStringVector* IFOweights; /**< list of relative detector-weights "w1, w2, w3, .." */

  REAL8 Freq;		/**< target-frequency */
  REAL8 dFreq;		/**< target-frequency offset */

  REAL8 Alpha;		/**< skyposition Alpha: radians, equatorial coords. */
  REAL8 dAlpha;		/**< skyposition Alpha offset */

  REAL8 Delta;		/**< skyposition Delta: radians, equatorial coords. */
  REAL8 dDelta;		/**< skyposition Delta offset */

  REAL8 f1dot;		/**< target 1. spindown-value df/dt */
  REAL8 df1dot;		/**< spindown df/dt offset */

  CHAR *ephemDir;	/**< directory of ephemeris-files */
  CHAR *ephemYear;	/**< year-range of ephemeris-file to use */

  REAL8 startTime;	/**< GPS start time of observation */
  REAL8 refTime;	/**< reference-time for spin-parameters fkdot */
  REAL8 duration;	/**< length of observation in seconds */
  INT4 numSteps;	/**< how many timesteps to use in Gauss-Legendre integration */

  REAL8 cosi;		/**< cos(iota) */
  REAL8 psi;		/**< polarization-angle psi */

  BOOLEAN printMotion;	/**< output orbital motion? */
  CHAR* outputMetric;	/**< filename to write metrics into */

  INT4 metricType;	/**< metric to compute: F-, phase-, orbital-, Ptole- metric */
  INT4 unitsType;	/**< which units to use for 't' and 'rX/c': SI vs natural */

  INT4 projection;     /**< project metric onto surface */
} UserVariables_t;


/*---------- empty structs for initializations ----------*/
ConfigVariables empty_ConfigVariables;
PulsarTimesParamStruc empty_PulsarTimesParamStruc;
UserVariables_t empty_UserVariables;
/* ---------- global variables ----------*/

extern int vrbflg;


/* ---------- local prototypes ---------- */
void initUserVars (LALStatus *status, UserVariables_t *uvar);
void InitCode (LALStatus *status, ConfigVariables *cfg, const UserVariables_t *uvar);
void getMultiPhaseDerivs (LALStatus *, MultiPhaseDerivs **derivs,
			  const MultiDetectorStateSeries *detStates,
			  const DopplerPoint *dopplerPoint,
			  PhaseType_t type,
			  const DopplerPoint *offsetUnits );

void InitEphemeris (LALStatus *, EphemerisData *edat, const CHAR *ephemDir, const CHAR *ephemYear, LIGOTimeGPS epoch, BOOLEAN isLISA);

int computeFstatMetric ( gsl_matrix *gF_ij, gsl_matrix *gFav_ij,
			 gsl_matrix *m1_ij, gsl_matrix *m2_ij, gsl_matrix *m3_ij,
			 ConfigVariables *cfg );
int computePhaseMetric ( gsl_matrix *g_ij, const PhaseDerivs *dphi, const REAL8Vector *GLweights );

int project_metric( gsl_matrix *ret_ij, gsl_matrix *g_ij, const UINT4 coordinate );
int outer_product ( gsl_matrix *ret_ij, const gsl_vector *u_i, const gsl_vector *v_j );
int symmetrize ( gsl_matrix *mat );
REAL8 quad_form ( const gsl_matrix *mat, const gsl_vector *vec );


void getPtolePosVel( PosVel_t *posvel, REAL8 tGPS, REAL8 tAutumn );

void XLALDestroyMultiPhaseDerivs ( MultiPhaseDerivs *mdPhi );

void FreeMem ( LALStatus *status, ConfigVariables *cfg );

void gauleg(double x1, double x2, double x[], double w[], int n);

/*============================================================
 * FUNCTION definitions
 *============================================================*/

int
main(int argc, char *argv[])
{
  LALStatus status = blank_status;
  ConfigVariables config = empty_ConfigVariables;
  gsl_matrix *m1_ij, *m2_ij, *m3_ij;
  gsl_matrix *gF_ij, *gFav_ij;
  REAL8 m1, m2, m3;
  REAL8 mF, mFav, disc, mMin, mMax;
  gsl_matrix *g_ij, *gFlat_ij, *gamma_ij;
  REAL8 mm;
  UserVariables_t uvar = empty_UserVariables;

  FILE *fpMetric = 0;
  PhaseType_t phaseType;
  MetricType_t startMetricType, stopMetricType, metricType;
  CHAR dummy[512];
  sprintf (dummy, "%s", lalGitID );
  sprintf (dummy, "%s", lalappsGitID );

  lalDebugLevel = 0;
  vrbflg = 1;	/* verbose error-messages */

  /* set LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;

  /* register user-variables */
  LAL_CALL (LALGetDebugLevel (&status, argc, argv, 'v'), &status);

  /* set log-level */
  LogSetLevel ( lalDebugLevel );

  LAL_CALL (initUserVars (&status, &uvar), &status);

  /* read cmdline & cfgfile  */
  LAL_CALL (LALUserVarReadAllInput (&status, argc,argv), &status);

  if (uvar.help) 	/* help requested: we're done */
    exit (0);

  if ( uvar.outputMetric )
    if ( (fpMetric = fopen ( uvar.outputMetric, "wb" )) == NULL )
      return FSTATMETRIC_EFILE;

  /* basic setup and initializations */
  gF_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  gFav_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  m1_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  m2_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  m3_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  g_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  gFlat_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  gamma_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );

  LAL_CALL ( InitCode(&status, &config, &uvar), &status);

  if ( uvar.metricType == METRIC_ALL )
    {
      startMetricType = METRIC_FSTAT;
      stopMetricType = METRIC_LAST;
    }
  else
    {
      startMetricType = uvar.metricType;
      stopMetricType = uvar.metricType + 1;
    }
  for ( metricType = startMetricType; metricType < stopMetricType; metricType ++ )
    {
      /* ----- compute phase derivatives ----- */
      if ( (metricType == METRIC_FSTAT ) || ( metricType == METRIC_PHASE ) )
	phaseType = PHASE_FULL;
      else if ( metricType == METRIC_ORBITAL )
	phaseType = PHASE_ORBITAL;
      else if ( metricType == METRIC_PTOLE )
	phaseType = PHASE_PTOLE;
      else
	phaseType = PHASE_NONE;
      if ( phaseType > PHASE_NONE )
	LAL_CALL ( getMultiPhaseDerivs (&status, &config.multidPhi, config.multiDetStates,
					&(config.dopplerPoint), phaseType, &(config.offsetUnits)), &status );

      switch ( metricType )
	{
	case METRIC_FSTAT:

	  if ( computeFstatMetric ( gF_ij, gFav_ij, m1_ij, m2_ij, m3_ij, &config ) )
	    {
	      printf ("\nSomething failed in computeFstatMetric() \n\n");
	      return -1;
	    }

          if ( uvar.projection > 0 )
            {
              project_metric ( gamma_ij, gF_ij, (uvar.projection - 1) );
              gsl_matrix_memcpy ( gF_ij, gamma_ij );
              project_metric ( gamma_ij, gFav_ij, (uvar.projection - 1) );
              gsl_matrix_memcpy ( gFav_ij, gamma_ij );
            }

	  mF   = quad_form ( gF_ij,   config.dopplerOffset );
	  mFav = quad_form ( gFav_ij, config.dopplerOffset );

	  m1 = quad_form ( m1_ij, config.dopplerOffset );
	  m2 = quad_form ( m2_ij, config.dopplerOffset );
	  m3 = quad_form ( m3_ij, config.dopplerOffset );

	  disc = mFav * mFav - ( m1 * m2 - m3 * m3 ) / config.Dd;
	  if ( disc < 0 )
	    {
	      CHAR *logstr = NULL;
	      LAL_CALL( LALUserVarGetLog (&status, &logstr, UVAR_LOGFMT_CMDLINE ), &status );
	      fprintf (stderr, "\nNegative discriminant %g, fudging this to plus...\n", disc);
	      fprintf (stderr, "Commandline was: %s\n\n", logstr );
	      LALFree ( logstr );
	      disc = - disc;
	    }
	  disc = sqrt(disc);

	  mMin = mFav - disc;
	  mMax = mFav + disc;

	  if ( fpMetric )
	    {
	      fprintf ( fpMetric, "\nA = %.16g;\nB = %.16g;\nC = %.16g;\nD = %.16g;\n",
			config.Ad, config.Bd, config.Cd, config.Dd );

	      fprintf ( fpMetric, "\ngF_ij = \\\n" ); XLALfprintfGSLmatrix ( fpMetric, METRIC_FORMAT,  gF_ij );
	      fprintf ( fpMetric, "\ngFav_ij = \\\n");XLALfprintfGSLmatrix ( fpMetric, METRIC_FORMAT,  gFav_ij );
	      fprintf ( fpMetric, "\nm1_ij = \\\n");  XLALfprintfGSLmatrix ( fpMetric, METRIC_FORMAT,  m1_ij );
	      fprintf ( fpMetric, "\nm2_ij = \\\n");  XLALfprintfGSLmatrix ( fpMetric, METRIC_FORMAT,  m2_ij );
	      fprintf ( fpMetric, "\nm3_ij = \\\n");  XLALfprintfGSLmatrix ( fpMetric, METRIC_FORMAT,  m3_ij );

	      fprintf ( fpMetric, "\nmF = %.16g;\nmFav = %.16g;\nmMin = %.16g;\nmMax = %.16g;\n\n",
			mF, mFav, mMin, mMax );
	    } /* if fpMetric */

	  break;

	  /* variants of the phase-metric */
	case METRIC_PHASE:
	case METRIC_ORBITAL:
	case METRIC_PTOLE:

	  /* only single-detector case allowed  */
	  if ( computePhaseMetric ( g_ij, config.multidPhi->data[0], config.GLweights) )
	    {
	      printf ("\nSomething failed in computePhaseMetric() \n\n");
	      return -1;
	    }

          if ( uvar.projection > 0 )
            {
              project_metric ( gamma_ij, g_ij, (uvar.projection - 1) );
              gsl_matrix_memcpy ( g_ij, gamma_ij );
            }

	  mm = quad_form ( g_ij, config.dopplerOffset );

	  if ( fpMetric )
	    {
	      const CHAR *gprefix, *mprefix;
	      if ( metricType == METRIC_PHASE ) {
		gprefix = "gPh_ij = \\\n"; mprefix = "mPh = ";
	      } else if ( metricType == METRIC_ORBITAL ) {
		gprefix = "gOrb_ij = \\\n"; mprefix = "mOrb = ";
	      } else if ( metricType == METRIC_PTOLE ) {
		gprefix = "gPtole_ij = \\\n"; mprefix = "mPtole = ";
	      }

              fprintf ( fpMetric, gprefix ); XLALfprintfGSLmatrix ( fpMetric, METRIC_FORMAT, g_ij );
              fprintf ( fpMetric, "\n%s %.16g;\n\n", mprefix, mm );

	    } /* if fpMetric */

	  break;

	case METRIC_FLAT:
	  {
	    REAL8 dntX, dntY, dom0, dom1; 	/* 'canonical' Doppler-variables */
	    gsl_vector *dopplerOffsetCanon;
	    REAL8 n1x, n1y, n1z, n2x, n2y, n2z, dnx, dny, dnz, dnX, dnY;
	    REAL8 sind1, sina1, cosd1, cosa1, sind2, sina2, cosd2, cosa2;
	    REAL8 sineps, coseps;
	    REAL8 Tspan = uvar.duration;
	    REAL8 Rorb = LAL_AU_SI;

	    /* ----- translate Doppler-offsets into 'canonical coords ----- */
	    if ( config.edat->leap < 0 )	/* signals that we're dealing with LISA */
	      {
		sineps = 0; coseps = 1;	/* already working in ecliptic coords */
	      }
	    else
	      {
		sineps = sin ( LAL_IEARTH ); coseps = cos ( LAL_IEARTH );
	      }

	    sind1 = sin(uvar.Delta); cosd1 = cos(uvar.Delta);
	    sina1 = sin(uvar.Alpha); cosa1 = cos(uvar.Alpha);

	    sind2 = sin(uvar.Delta + uvar.dDelta); cosd2 = cos(uvar.Delta + uvar.dDelta);
	    sina2 = sin(uvar.Alpha + uvar.dAlpha); cosa2 = cos(uvar.Alpha + uvar.dAlpha);


	    n1x = cosa1 * cosd1;
	    n1y = sina1 * cosd1;
	    n1z = sind1;

	    n2x = cosa2 * cosd2;
	    n2y = sina2 * cosd2;
	    n2z = sind2;
	    /* sky-pos offset in input sky-coords */
	    dnx = n2x - n1x;
	    dny = n2y - n1y;
	    dnz = n2z - n1z;

	    /* translate into ecliptic corrds */
	    dnX = dnx;
	    dnY = coseps * dny + sineps * dnz;

	    /* sky-pos offset in canonical units */
	    dntX = LAL_TWOPI * uvar.Freq * Rorb / LAL_C_SI * dnX;
	    dntY = LAL_TWOPI * uvar.Freq * Rorb / LAL_C_SI * dnY;

	    /* spin-offsets in canonical units */
	    dom0 = 1.0 * LAL_TWOPI * Tspan * uvar.dFreq;	/* / (s+1)! */
	    dom1 = 0.5 * LAL_TWOPI * Tspan * Tspan * uvar.df1dot;

	    dopplerOffsetCanon = gsl_vector_calloc ( METRIC_DIM );
	    gsl_vector_set ( dopplerOffsetCanon, 0, dom0 );
	    gsl_vector_set ( dopplerOffsetCanon, 1, dntX );
	    gsl_vector_set ( dopplerOffsetCanon, 2, dntY );
	    gsl_vector_set ( dopplerOffsetCanon, 3, dom1 );

	    LogPrintf (LOG_DETAIL, "Offsets: dom0 = %.15g, dntX = %.15g, dntY = %.15g, dom1 = %.15g\n",
		       dom0, dntX, dntY, dom1 );

	    if ( 0 != XLALFlatMetricCW ( gFlat_ij, config.refTime, config.startTime, Tspan, config.edat ) )
	      {
		LogPrintf ( LOG_CRITICAL, "XLALFlatMetricCW() failed!\n");
		return -1;
	      }

	    mm = quad_form ( gFlat_ij, dopplerOffsetCanon );
	    if ( fpMetric )
	      {
		const CHAR *gprefix = "gFlat_ij = \\\n";
		const CHAR *mprefix = "mFlat = ";

		fprintf ( fpMetric, gprefix); XLALfprintfGSLmatrix ( fpMetric, METRIC_FORMAT, gFlat_ij );
		fprintf ( fpMetric, "\n%s %.16g;\n\n", mprefix, mm );
	      } /* if fpMetric */

	  }
	  break;

	default:
	  LALPrintError("\nInvalid metric-number '%d'.\n\n", metricType );
	  return -1;
	  break;
	} /* switch ( metricType ) */

      XLALDestroyMultiPhaseDerivs ( config.multidPhi );
      config.multidPhi = NULL;

    } /* for metricType < lastMetricType */

  if ( fpMetric )
    fclose(fpMetric);

  /* ----- done: free all memory */
  gsl_matrix_free ( gF_ij );
  gsl_matrix_free ( gFav_ij );
  gsl_matrix_free ( g_ij );
  gsl_matrix_free ( gFlat_ij );
  gsl_matrix_free ( gamma_ij );

  gsl_matrix_free ( m1_ij );
  gsl_matrix_free ( m2_ij );
  gsl_matrix_free ( m3_ij );

  LAL_CALL (FreeMem(&status, &config), &status);

  LALCheckMemoryLeaks();

  return 0;
} /* main */


/* calculate Fstat-metric components m1_ij, m2_ij, m3_ij, which must
 * be allocated 4x4 matrices.
 *
 * Return 0 = OK, -1 on error.
 */
int
computeFstatMetric ( gsl_matrix *gF_ij, gsl_matrix *gFav_ij,
		     gsl_matrix *m1_ij, gsl_matrix *m2_ij, gsl_matrix *m3_ij,
		     ConfigVariables *cfg )
{
  UINT4 i;
  REAL8 A, B, C, D;
  UINT4 X, numDet;
  REAL8 AMA;

  gsl_matrix *P1_ij, *P2_ij, *P3_ij;
  gsl_vector *a2_dPhi_i, *b2_dPhi_i, *ab_dPhi_i;

  gsl_matrix *Q1_ij, *Q2_ij, *Q3_ij;

  gsl_vector *dPhi_i;
  gsl_matrix *mat1, *mat2;
  gsl_vector *vec;

  gsl_matrix *a2_a2, *a2_b2, *a2_ab;
  gsl_matrix *b2_b2, *b2_ab;
  gsl_matrix *ab_ab;

  numDet = cfg->multiDetStates->length;

  /* ----- check input ----- */
  if ( !m1_ij || !m2_ij || !m3_ij )
    return -1;

  if ( (m1_ij->size1 != METRIC_DIM) || (m1_ij->size2 != METRIC_DIM) )
    return -1;
  if ( (m2_ij->size1 != METRIC_DIM) || (m2_ij->size2 != METRIC_DIM) )
    return -1;
  if ( (m3_ij->size1 != METRIC_DIM) || (m3_ij->size2 != METRIC_DIM) )
    return -1;

  if ( cfg->multiAMcoe->length != numDet )
    return -1;
  if ( cfg->multidPhi->length != numDet )
    return -1;

  /* ----- allocate matrices/vectors ----- */
  dPhi_i = gsl_vector_calloc ( METRIC_DIM );

  mat1 = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  mat2 = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  vec = gsl_vector_calloc (  METRIC_DIM );

  P1_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  P2_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  P3_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );

  a2_dPhi_i = gsl_vector_calloc ( METRIC_DIM );
  b2_dPhi_i = gsl_vector_calloc ( METRIC_DIM );
  ab_dPhi_i = gsl_vector_calloc ( METRIC_DIM );

  Q1_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  Q2_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  Q3_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );

  a2_a2 = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  a2_b2 = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  a2_ab = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );

  b2_b2 = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  b2_ab = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );

  ab_ab = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );

  /* ----- calculate averages ----- */
  for ( X=0; X < numDet; X ++ )
    {
      UINT4 numSteps = cfg->multiDetStates->data[X]->length;
      AMCoeffs *amcoe = cfg->multiAMcoe->data[X];
      PhaseDerivs *dPhi = cfg->multidPhi->data[X];
      gsl_matrix *P1_Xij, *P2_Xij, *P3_Xij;
      gsl_vector *a2_dPhi_Xi, *b2_dPhi_Xi, *ab_dPhi_Xi;

      P1_Xij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
      P2_Xij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
      P3_Xij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );

      a2_dPhi_Xi = gsl_vector_calloc ( METRIC_DIM );
      b2_dPhi_Xi = gsl_vector_calloc ( METRIC_DIM );
      ab_dPhi_Xi = gsl_vector_calloc ( METRIC_DIM );

      for ( i = 0; i < numSteps; i ++ )
	{
	  REAL8 a = amcoe->a->data[i];
	  REAL8 b = amcoe->b->data[i];
	  REAL8 a2 = SQ(a);
	  REAL8 b2 = SQ(b);
	  REAL8 ab = a * b;

	  gsl_vector_set ( dPhi_i, 0, dPhi->dFreq->data[i] );
	  gsl_vector_set ( dPhi_i, 1, dPhi->dAlpha->data[i]);
	  gsl_vector_set ( dPhi_i, 2, dPhi->dDelta->data[i]);
	  gsl_vector_set ( dPhi_i, 3, dPhi->df1dot->data[i]);

	  outer_product(mat1, dPhi_i, dPhi_i );

	  /* ----- P1_ij ----- */
	  gsl_matrix_memcpy (mat2, mat1);
	  gsl_matrix_scale ( mat2, a2 );
	  gsl_matrix_add ( P1_Xij, mat2 );

	  /* ----- P2_ij ----- */
	  gsl_matrix_memcpy (mat2, mat1);
	  gsl_matrix_scale ( mat2, b2 );
	  gsl_matrix_add ( P2_Xij, mat2 );

	  /* ----- P3_ij ----- */
	  gsl_matrix_memcpy (mat2, mat1);
	  gsl_matrix_scale ( mat2, ab );
	  gsl_matrix_add ( P3_Xij, mat2 );

	  /* ----- a2_dPhi_i ----- */
	  gsl_vector_memcpy ( vec, dPhi_i );
	  gsl_vector_scale ( vec, a2 );
	  gsl_vector_add ( a2_dPhi_Xi, vec );

	  /* ----- b2_dPhi_i ----- */
	  gsl_vector_memcpy ( vec, dPhi_i );
	  gsl_vector_scale ( vec, b2 );
	  gsl_vector_add ( b2_dPhi_Xi, vec );

	  /* ----- ab_dPhi_i ----- */
	  gsl_vector_memcpy ( vec, dPhi_i );
	  gsl_vector_scale ( vec, ab );
	  gsl_vector_add ( ab_dPhi_Xi, vec );

	} /* for i < numSteps */

      gsl_matrix_add ( P1_ij, P1_Xij );
      gsl_matrix_add ( P2_ij, P2_Xij );
      gsl_matrix_add ( P3_ij, P3_Xij );

      gsl_vector_add ( a2_dPhi_i, a2_dPhi_Xi );

      gsl_vector_add ( b2_dPhi_i, b2_dPhi_Xi );

      gsl_vector_add ( ab_dPhi_i, ab_dPhi_Xi );

      gsl_matrix_free ( P1_Xij );
      gsl_matrix_free ( P2_Xij );
      gsl_matrix_free ( P3_Xij );

      gsl_vector_free ( a2_dPhi_Xi );
      gsl_vector_free ( b2_dPhi_Xi );
      gsl_vector_free ( ab_dPhi_Xi );

    } /* for X < numDet */

  /* ---------- composite quantities ---------- */
  A = cfg->multiAMcoe->Mmunu.Ad ;
  B = cfg->multiAMcoe->Mmunu.Bd ;
  C = cfg->multiAMcoe->Mmunu.Cd ;
  D = A * B - C*C;
  cfg->Ad = A;
  cfg->Bd = B;
  cfg->Cd = C;
  cfg->Dd = D;

  outer_product (a2_a2, a2_dPhi_i, a2_dPhi_i );
  outer_product (a2_b2, a2_dPhi_i, b2_dPhi_i );
  outer_product (a2_ab, a2_dPhi_i, ab_dPhi_i );

  outer_product (b2_b2, b2_dPhi_i, b2_dPhi_i );
  outer_product (b2_ab, b2_dPhi_i, ab_dPhi_i );

  outer_product (ab_ab, ab_dPhi_i, ab_dPhi_i );

  /* ----- Q1_ij ----- */
  gsl_matrix_memcpy ( mat1, ab_ab );
  gsl_matrix_scale ( mat1, A/D );
  gsl_matrix_memcpy ( Q1_ij, mat1 );	/*  = (A/D)<ab_dPhi_i><ab_dPhi_j> */

  gsl_matrix_memcpy ( mat1, a2_a2 );
  gsl_matrix_scale ( mat1, B/D );
  gsl_matrix_add ( Q1_ij, mat1 );	/*  + (B/D)<a2_dPhi_i><a2_dPhi_j> */

  gsl_matrix_memcpy ( mat1, a2_ab );
  gsl_matrix_scale ( mat1, - 2.0 * C/D );
  gsl_matrix_add ( Q1_ij, mat1 );	/*  -2(C/D)<a2_dPhi_i><ab_dPhi_j> */

  symmetrize ( Q1_ij );		/* (i,j) */

  /* ----- Q2_ij ----- */
  gsl_matrix_memcpy ( mat1, b2_b2 );
  gsl_matrix_scale ( mat1, A/D );
  gsl_matrix_memcpy ( Q2_ij, mat1 );	/*  = (A/D)<b2_dPhi_i><b2_dPhi_j> */

  gsl_matrix_memcpy ( mat1, ab_ab );
  gsl_matrix_scale ( mat1, B/D );
  gsl_matrix_add ( Q2_ij, mat1 );	/*  + (B/D)<ab_dPhi_i><ab_dPhi_j> */

  gsl_matrix_memcpy ( mat1, b2_ab );
  gsl_matrix_scale ( mat1, - 2.0 * C/D );
  gsl_matrix_add ( Q2_ij, mat1 );	/*  -2(C/D)<b2_dPhi_i><ab_dPhi_j> */

  symmetrize ( Q2_ij );		/* (i,j) */

  /* ----- Q3_ij ----- */
  gsl_matrix_memcpy ( mat1, b2_ab );
  gsl_matrix_scale ( mat1, A/D );
  gsl_matrix_memcpy ( Q3_ij, mat1 );	/*  = (A/D)<b2_dPhi_i><ab_dPhi_j> */

  gsl_matrix_memcpy ( mat1, a2_ab );
  gsl_matrix_scale ( mat1, B/D );
  gsl_matrix_add ( Q3_ij, mat1 );	/*  + (B/D)<a2_dPhi_i><ab_dPhi_j> */

  gsl_matrix_memcpy ( mat1, a2_b2 );
  gsl_matrix_add ( mat1, ab_ab );
  gsl_matrix_scale ( mat1, - 1.0 * C/D );
  gsl_matrix_add ( Q3_ij, mat1 );	/*  -(C/D)(<a2_dPhi_i><b2_dPhi_j> + <ab_dPhi_i><ab_dPhi_j>)*/

  symmetrize ( Q3_ij );		/* (i,j) */

  /* ===== final matrics: m1_ij, m2_ij, m3_ij ===== */
  gsl_matrix_memcpy ( m1_ij, P1_ij );
  gsl_matrix_sub ( m1_ij, Q1_ij );

  gsl_matrix_memcpy ( m2_ij, P2_ij );
  gsl_matrix_sub ( m2_ij, Q2_ij );

  gsl_matrix_memcpy ( m3_ij, P3_ij );
  gsl_matrix_sub ( m3_ij, Q3_ij );


  /* ===== full F-metric gF_ij */
  AMA = A  * cfg->Al1 + B  * cfg->Al2 + 2.0 * C  * cfg->Al3;

  gsl_matrix_memcpy (gF_ij, m1_ij );
  gsl_matrix_scale ( gF_ij, cfg->Al1 );	/* alpha1 m^1_ij */

  gsl_matrix_memcpy (mat1, m2_ij );
  gsl_matrix_scale ( mat1, cfg->Al2 );
  gsl_matrix_add ( gF_ij, mat1 );	/* + alpha2 m^2_ij */

  gsl_matrix_memcpy (mat1, m3_ij );
  gsl_matrix_scale ( mat1, 2.0 * cfg->Al3 );
  gsl_matrix_add ( gF_ij, mat1 );	/* + 2 * alpha3 m^3_ij */

  gsl_matrix_scale ( gF_ij, 1.0 / AMA );

  /* ===== averaged F-metric gFav_ij */
  gsl_matrix_memcpy (gFav_ij, m1_ij );
  gsl_matrix_scale ( gFav_ij, B );	/* B m^1_ij */

  gsl_matrix_memcpy (mat1, m2_ij );
  gsl_matrix_scale ( mat1, A );
  gsl_matrix_add ( gFav_ij, mat1 );	/* + A m^2_ij */

  gsl_matrix_memcpy (mat1, m3_ij );
  gsl_matrix_scale ( mat1, - 2.0 * C );
  gsl_matrix_add ( gFav_ij, mat1 );	/* - 2C m^3_ij */

  gsl_matrix_scale ( gFav_ij, 1.0 / ( 2.0 * D ) ); /* 1/ (2D) */


  /* ----- free memory ----- */
  gsl_vector_free ( dPhi_i );

  gsl_matrix_free ( mat1 );
  gsl_matrix_free ( mat2 );
  gsl_vector_free ( vec );

  gsl_matrix_free ( P1_ij );
  gsl_matrix_free ( P2_ij );
  gsl_matrix_free ( P3_ij );

  gsl_vector_free ( a2_dPhi_i );
  gsl_vector_free ( b2_dPhi_i );
  gsl_vector_free ( ab_dPhi_i );

  gsl_matrix_free ( Q1_ij );
  gsl_matrix_free ( Q2_ij );
  gsl_matrix_free ( Q3_ij );

  gsl_matrix_free ( a2_a2 );
  gsl_matrix_free ( a2_b2 );
  gsl_matrix_free ( a2_ab );

  gsl_matrix_free ( b2_b2 );
  gsl_matrix_free ( b2_ab );

  gsl_matrix_free ( ab_ab );


  return 0;

} /* computeFstatMetric() */

/* calculate pure Phase-metric gij, which must be allocated 4x4 matrix.
 *
 * Return 0 = OK, -1 on error.
 */
int
computePhaseMetric ( gsl_matrix *g_ij, const PhaseDerivs *dPhi, const REAL8Vector *GLweights )
{
  UINT4 i;
  UINT4 numSteps = GLweights->length;

  gsl_vector *dPhi_i;
  gsl_matrix *dPhi_i_dPhi_j;

  gsl_matrix *aPhi_ij;
  gsl_vector *aPhi_i;
  gsl_matrix *aPhi_i_aPhi_j;

  /* ----- check input ----- */
  if ( !g_ij  )
    return -1;

  if ( (g_ij->size1 != METRIC_DIM) || (g_ij->size2 != METRIC_DIM) )
    return -1;

  if ( dPhi->dAlpha->length != numSteps )
    return -1;

  /* ----- allocate matrices/vectors ----- */
  dPhi_i = gsl_vector_calloc ( METRIC_DIM );
  dPhi_i_dPhi_j = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );

  aPhi_i = gsl_vector_calloc ( METRIC_DIM );
  aPhi_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  aPhi_i_aPhi_j = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );

  /* ----- calculate averages using Gauss-Legendre integration ----- */
  for ( i = 0; i < numSteps; i ++ )
    {
      REAL8 wi = GLweights->data[i];
      gsl_vector_set ( dPhi_i, 0, dPhi->dFreq->data[i] );
      gsl_vector_set ( dPhi_i, 1, dPhi->dAlpha->data[i]);
      gsl_vector_set ( dPhi_i, 2, dPhi->dDelta->data[i]);
      gsl_vector_set ( dPhi_i, 3, dPhi->df1dot->data[i]);

      outer_product(dPhi_i_dPhi_j, dPhi_i, dPhi_i );

      /* Gauss-Legendre integration: weighted sum */
      gsl_vector_scale ( dPhi_i, wi );
      gsl_matrix_scale ( dPhi_i_dPhi_j, wi );

      /* ----- av_dPhi_i ----- */
      gsl_vector_add ( aPhi_i, dPhi_i );

      /* ----- av_dPhi_i_dPhi_j ----- */
      gsl_matrix_add ( aPhi_ij, dPhi_i_dPhi_j );

    } /* for i < numSteps */

  /* ---------- composite quantities ---------- */
  outer_product (aPhi_i_aPhi_j, aPhi_i, aPhi_i );

  /* ===== final answer: m1_ij, m2_ij, m3_ij ===== */
  gsl_matrix_memcpy ( g_ij, aPhi_ij );
  gsl_matrix_sub ( g_ij, aPhi_i_aPhi_j );

  /* ----- free memory ----- */
  gsl_vector_free ( dPhi_i );
  gsl_matrix_free ( dPhi_i_dPhi_j );
  gsl_matrix_free ( aPhi_ij );
  gsl_vector_free ( aPhi_i );
  gsl_matrix_free ( aPhi_i_aPhi_j );

  return 0;

} /* computePhaseMetric() */


/** register all "user-variables" */
void
initUserVars (LALStatus *status, UserVariables_t *uvar)
{
  INITSTATUS( status, "initUserVars", rcsid );
  ATTATCHSTATUSPTR (status);

  /* set a few defaults */
  uvar->help = FALSE;

#define EPHEM_YEARS  "00-04"
  uvar->ephemYear = (CHAR*)LALCalloc (1, strlen(EPHEM_YEARS)+1);
  strcpy (uvar->ephemYear, EPHEM_YEARS);

#define DEFAULT_EPHEMDIR "env LAL_DATA_PATH"
  uvar->ephemDir = (CHAR*)LALCalloc (1, strlen(DEFAULT_EPHEMDIR)+1);
  strcpy (uvar->ephemDir, DEFAULT_EPHEMDIR);

  uvar->Freq = 100;
  uvar->f1dot = 0.0;

  uvar->startTime = 714180733;
  uvar->refTime   = uvar->startTime;
  uvar->duration = 10 * 3600;
  uvar->numSteps = 2000;

  uvar->printMotion = FALSE;

  uvar->metricType = METRIC_ALL;
  uvar->unitsType = UNITS_SI;

  uvar->projection = 0;

  /* register all our user-variable */

  LALregBOOLUserStruct(status,	help,		'h', UVAR_HELP,		"Print this help/usage message");
  LALregLISTUserStruct(status,	IFOs,		'I', UVAR_REQUIRED, 	"Comma-separated list of detectors, eg. \"H1,H2,L1,G1, ...\" ");
  LALregLISTUserStruct(status,	IFOweights,	 0, UVAR_OPTIONAL, 	"Comma-separated list of relative noise-weights, eg. \"w1,w2,w3,..\" ");
  LALregREALUserStruct(status,	Alpha,		'a', UVAR_OPTIONAL,	"skyposition Alpha in radians, equatorial coords.");
  LALregREALUserStruct(status,	dAlpha,		 0, UVAR_OPTIONAL,	"skyposition offset in Alpha");
  LALregREALUserStruct(status,	Delta, 		'd', UVAR_OPTIONAL,	"skyposition Delta in radians, equatorial coords.");
  LALregREALUserStruct(status,	dDelta, 	 0, UVAR_OPTIONAL,	"skyposition offset in Delta");
  LALregREALUserStruct(status,	Freq, 		'f', UVAR_OPTIONAL, 	"target frequency");
  LALregREALUserStruct(status,	dFreq, 		 0, UVAR_OPTIONAL, 	"target frequency offset");
  LALregREALUserStruct(status,	f1dot, 		's', UVAR_OPTIONAL, 	"first spindown-value df/dt");
  LALregREALUserStruct(status,	df1dot, 	 0, UVAR_OPTIONAL, 	"first spindown-value offset");
  LALregREALUserStruct(status, 	startTime,      't', UVAR_OPTIONAL, 	"GPS start time of observation");
  LALregREALUserStruct(status, 	refTime,      	 0, UVAR_OPTIONAL, 	"reference time for spin-parameters [Default = startTime]");
  LALregREALUserStruct(status,  duration,	'T', UVAR_OPTIONAL,	"Alternative: Duration of observation in seconds");
  LALregINTUserStruct(status,   numSteps,	 0,  UVAR_OPTIONAL,	"Order of Gauss-Legendre quadrature to use");
  LALregSTRINGUserStruct(status, ephemDir,       'E', UVAR_OPTIONAL, 	"Directory where Ephemeris files are located");
  LALregSTRINGUserStruct(status, ephemYear,      'y', UVAR_OPTIONAL, 	"Year (or range of years) of ephemeris files to be used");
  LALregREALUserStruct(status, 	cosi,	 	 0, UVAR_OPTIONAL,	"Pulsar orientation-angle cos(iota) [-1,1]" );
  LALregREALUserStruct(status,	psi,		 0, UVAR_OPTIONAL,	"Wave polarization-angle psi [0, pi]" );
  LALregBOOLUserStruct(status,	printMotion,	 0, UVAR_OPTIONAL,	"Output the orbital motion for integration-steps.");
  LALregSTRINGUserStruct(status, outputMetric,	 0, UVAR_OPTIONAL,	"Output the metric components (in octave format) into this file.");
  LALregINTUserStruct(status,  	metricType,	 0,  UVAR_OPTIONAL,	"Type of metric: 0=ALL, 1=mF, 2=mPh, 3=mOrb, 4=mPtole, 5=flat");
  LALregINTUserStruct(status,  	unitsType,	 0,  UVAR_OPTIONAL,	"Which units to use for metric components: 0=SI, 1=natural (order-1)");
  LALregINTUserStruct(status,   projection,      0,  UVAR_OPTIONAL,     "Coordinate of metric projection: 0=none, 1=f, 2=Alpha, 3=Delta, 4=f1dot");

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* initUserVars() */


/** basic initializations: set-up 'ConfigVariables'
 */
void
InitCode (LALStatus *status, ConfigVariables *cfg, const UserVariables_t *uvar)
{
  LIGOTimeGPSVector *GLtimestamps = NULL;
  REAL8Vector *detWeights;
  UINT4 X, numDet;

  INITSTATUS (status, "InitCode", rcsid);
  ATTATCHSTATUSPTR (status);

  /* ----- determine start-time from user-input */
  XLALGPSSetREAL8( &(cfg->startTime), uvar->startTime );
  /* if refTime wasn't set, use startTime ! */
  if ( LALUserVarWasSet ( &uvar->refTime ) )
    XLALGPSSetREAL8( &(cfg->refTime), uvar->refTime );
  else
    XLALGPSSetREAL8( &(cfg->refTime), uvar->startTime );


  { /* ----- load ephemeris-data ----- */
    CHAR *ephemDir;
    BOOLEAN isLISA = FALSE;

    cfg->edat = LALCalloc(1, sizeof(EphemerisData));
    if ( LALUserVarWasSet ( &uvar->ephemDir ) )
      ephemDir = uvar->ephemDir;
    else
      ephemDir = NULL;

    /* hack: if first detector is LISA, we load MLDC-ephemeris instead of 'earth' files */
    if ( uvar->IFOs->data[0][0] == 'Z' )
      isLISA = TRUE;

    TRY( InitEphemeris (status->statusPtr, cfg->edat, ephemDir, uvar->ephemYear, cfg->startTime, isLISA ), status);
  }

  /* ----- which units to use? 'natural' or SI ----- */
  if ( (cfg->offsetUnits.fkdot = XLALCreateREAL8Vector ( NUM_SPINS )) == NULL ) {
    ABORT (status, FSTATMETRIC_EMEM, FSTATMETRIC_MSGEMEM);
  }
  cfg->offsetUnits.skypos.system = COORDINATESYSTEM_EQUATORIAL;
  if ( uvar->unitsType == UNITS_NATURAL )
    {
      cfg->offsetUnits.fkdot->data[0] = 1.0 / ( LAL_TWOPI * uvar->duration );
      cfg->offsetUnits.fkdot->data[1] = 2.0 / ( LAL_TWOPI * SQ( uvar->duration ) );
      cfg->offsetUnits.skypos.longitude = 1.0;
      cfg->offsetUnits.skypos.latitude = 1.0;
    }
  else /* SI-units */
    {
      cfg->offsetUnits.fkdot->data[0] = 1.0;
      cfg->offsetUnits.fkdot->data[1] = 1.0;
      cfg->offsetUnits.skypos.longitude = 1.0;
      cfg->offsetUnits.skypos.latitude = 1.0;
    }

  /* ----- get parameter-space point from user-input) */
  {
    /* NOTE: internally we always deal with a fixed number of spins */
    if ( (cfg->dopplerPoint.fkdot = XLALCreateREAL8Vector ( NUM_SPINS )) == NULL ) {
      ABORT (status, FSTATMETRIC_EMEM, FSTATMETRIC_MSGEMEM);
    }

    cfg->dopplerPoint.fkdot->data[0] = uvar->Freq;
    cfg->dopplerPoint.fkdot->data[1] = uvar->f1dot;

    cfg->dopplerPoint.skypos.system = COORDINATESYSTEM_EQUATORIAL;
    cfg->dopplerPoint.skypos.longitude = uvar->Alpha;
    cfg->dopplerPoint.skypos.latitude = uvar->Delta;
  } /* get parameter-space point */

  /* ----- set Doppler-offset (take care of metric units) ----- */
  {
    REAL8 dFreq, dAlpha, dDelta, df1dot;

    dFreq 	= uvar->dFreq   / cfg->offsetUnits.fkdot->data[0];
    dAlpha	= uvar->dAlpha  / cfg->offsetUnits.skypos.longitude;
    dDelta	= uvar->dDelta  / cfg->offsetUnits.skypos.latitude;
    df1dot	= uvar->df1dot  / cfg->offsetUnits.fkdot->data[1];

    cfg->dopplerOffset = gsl_vector_calloc ( METRIC_DIM );
    gsl_vector_set ( cfg->dopplerOffset, 0, dFreq );
    gsl_vector_set ( cfg->dopplerOffset, 1, dAlpha );
    gsl_vector_set ( cfg->dopplerOffset, 2, dDelta );
    gsl_vector_set ( cfg->dopplerOffset, 3, df1dot );
  } /* set Doppler-offset */

  /* ----- construct Gauss-Legendre timestamps and corresponding weights
   * for computing the integrals by Gauss-Legendre quadrature
   */
  {
    UINT4 i;
    REAL8Vector *ti;	/* temporary */
    REAL8 Tinv = 1.0 / uvar->duration;

    if ( (cfg->GLweights = XLALCreateREAL8Vector ( uvar->numSteps )) == NULL ) {
      ABORT (status, FSTATMETRIC_EMEM, FSTATMETRIC_MSGEMEM );
    }
    if ( ( ti = XLALCreateREAL8Vector ( uvar->numSteps )) == NULL ) {
      ABORT (status, FSTATMETRIC_EMEM, FSTATMETRIC_MSGEMEM );
    }
    if ( (GLtimestamps = XLALCreateTimestampVector ( uvar->numSteps )) == NULL ) {
      ABORT (status, FSTATMETRIC_EMEM, FSTATMETRIC_MSGEMEM );
    }

    /* compute Gauss-Legendre roots, timestamps and associated weights */
    gauleg(uvar->startTime, uvar->startTime+uvar->duration, ti->data, cfg->GLweights->data, uvar->numSteps);

    /* convert REAL8-times into LIGOTimeGPS-times */
    for ( i=0; i < (UINT4)uvar->numSteps; i ++ )
      {
	XLALGPSSetREAL8 ( &(GLtimestamps->data[i]), ti->data[i] );
	cfg->GLweights->data[i] *= Tinv;
      }

    XLALDestroyREAL8Vector ( ti );

  } /* setup time-series of GL-timestamps */

  /* ----- initialize IFOs and (Multi-)DetectorStateSeries  ----- */
  {
    LALDetector *ifo = NULL;

    numDet = uvar->IFOs->length;

    if ( (cfg->multiDetStates = LALCalloc ( 1, sizeof( *(cfg->multiDetStates) ))) == NULL ) {
      ABORT (status, FSTATMETRIC_EMEM, FSTATMETRIC_MSGEMEM);
    }
    if ( (cfg->multiDetStates->data = LALCalloc (numDet, sizeof( *(cfg->multiDetStates->data))))==NULL){
      ABORT (status, FSTATMETRIC_EMEM, FSTATMETRIC_MSGEMEM);
    }
    cfg->multiDetStates->length = numDet;
    cfg->multiDetStates->Tspan = uvar->duration;

    for ( X=0; X < numDet; X ++ )
      {
	if ( ( ifo = XLALGetSiteInfo ( uvar->IFOs->data[X] ) ) == NULL ) {
	  LALPrintError("\nFailed to get site-info for IFO '%s'\n\n", uvar->IFOs->data[X] );
	  ABORT (status, FSTATMETRIC_EINPUT, FSTATMETRIC_MSGEINPUT);
	}
	/* obtain detector positions and velocities, together with LMSTs */
	TRY (LALGetDetectorStates(status->statusPtr, &(cfg->multiDetStates->data[X]), GLtimestamps,
				  ifo, cfg->edat, 0 ), status);
	LALFree ( ifo );

      } /* for i < numDet */

  } /* init multi-DetStateSeries */


  /* ----- get relative detector-weights from user ---------- */
  if ( (detWeights = XLALCreateREAL8Vector ( numDet )) == NULL ) {
    ABORT (status, FSTATMETRIC_EMEM, FSTATMETRIC_MSGEMEM);
  }
  for ( X=0; X < numDet ; X ++ )
    detWeights->data[X] = 1.0;	/* default: equal weights */

  if ( uvar->IFOweights )
    {
      if ( uvar->IFOweights->length != numDet )
	{
	  LALPrintError ("\nNumber of IFOweights must agree with IFOs if given!\n\n");
	  ABORT (status, FSTATMETRIC_EINPUT, FSTATMETRIC_MSGEINPUT);
	}
      for ( X=0; X < numDet ; X ++ )
	{
	  if ( 1 != sscanf ( uvar->IFOweights->data[X], "%lf", &(detWeights->data[X])) )
	    {
	      LALPrintError ("\nFailed to parse noise-weight '%s' into float.\n\n",
			     uvar->IFOweights->data[X] );
	      ABORT (status, FSTATMETRIC_EINPUT, FSTATMETRIC_MSGEINPUT);
	    }
	} /* for X < numDet */

    } /* if uvar->IFOweights */

  /* ---------- combine relative detector-weights with GL-weights ----------*/
  {
    MultiNoiseWeights *tmp;
    if ( (tmp = LALCalloc(1, sizeof(*tmp))) == NULL ){
      ABORT (status,  FSTATMETRIC_EMEM,  FSTATMETRIC_MSGEMEM);
    }
    tmp->length = numDet;

    if ( (tmp->data = LALCalloc( numDet, sizeof(*(tmp->data)))) == NULL) {
      ABORT (status,  FSTATMETRIC_EMEM,  FSTATMETRIC_MSGEMEM);
    }

    for ( X = 0; X < numDet; X ++)
      {
	UINT4 alpha;

	/* create k^th weights vector */
	if ( (tmp->data[X] = XLALCreateREAL8Vector ( uvar->numSteps )) == NULL ) {
	  ABORT (status,  FSTATMETRIC_EMEM,  FSTATMETRIC_MSGEMEM);
	}

	/* set all weights to detectorWeight * GLweight */
	for ( alpha = 0; alpha < (UINT4)uvar->numSteps; alpha ++ )
	  tmp->data[X]->data[alpha] = detWeights->data[X] * cfg->GLweights->data[alpha];

      } /* for X < numDet */

    cfg->multiNoiseWeights = tmp;

  } /* ----- construct multiNoiseWeights */

  TRY ( LALGetMultiAMCoeffs (status->statusPtr, &(cfg->multiAMcoe), cfg->multiDetStates,
			     cfg->dopplerPoint.skypos ), status );

  if ( XLALWeighMultiAMCoeffs( cfg->multiAMcoe, cfg->multiNoiseWeights ) != XLAL_SUCCESS )
    {
      LALPrintError ( "\nSomething failed in XLALWeighMultiAMCoeffs() ...\n\n");
      ABORT ( status, FSTATMETRIC_EINPUT, FSTATMETRIC_MSGEINPUT );
    }

  /* ----- compute amplitude-factors alpha1, alpha2, alpha3 ----- */
  {
    REAL8 Aplus = 0.5 * ( 1.0 + SQ(uvar->cosi) );
    REAL8 Across = uvar->cosi;
    REAL8 cos2psi = cos(2.0 * uvar->psi );
    REAL8 sin2psi = sin(2.0 * uvar->psi );
    cfg->Al1 = SQ(Aplus) * SQ( cos2psi ) + SQ(Across) * SQ(sin2psi);
    cfg->Al2 = SQ(Aplus) * SQ( sin2psi ) + SQ(Across) * SQ(cos2psi);
    cfg->Al3 = ( SQ(Aplus) - SQ(Across) ) * sin2psi * cos2psi ;
  }

  /* free temporary memory */
  XLALDestroyREAL8Vector ( detWeights );
  XLALDestroyTimestampVector ( GLtimestamps);

  DETATCHSTATUSPTR (status);
  RETURN (status);


} /* InitFStat() */


/** calculate the phase-derivatives \f$\partial_i \phi \f$ for the
 * time-series detStates and the given doppler-point.
 * Has the option of using only the orbital part of the phase (PHASE_ORBITAL)
 * or the full-phase (PHASE_FULL).
 *
 * returned PhaseDerivs is allocated in here.
 */
void
getMultiPhaseDerivs (LALStatus *status,
		     MultiPhaseDerivs **multidPhi,
		     const MultiDetectorStateSeries *multiDetStates,
		     const DopplerPoint *dopplerPoint,
		     PhaseType_t phaseType,
		     const DopplerPoint *offsetUnits )
{
  UINT4 numDet, X;
  UINT4 i;
  REAL8 vn[3];		/* unit-vector pointing to source in Cart. equatorial coord. */
  REAL8 Alpha, Delta;
  REAL8 refTime;
  LIGOTimeGPS refTimeGPS;
  PulsarTimesParamStruc times = empty_PulsarTimesParamStruc;
  MultiPhaseDerivs *mdPhi = NULL;
  REAL8 TspanInv;

  INITSTATUS (status, "getMultiPhaseDerivs", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT ( dopplerPoint, status, FSTATMETRIC_ENULL, FSTATMETRIC_MSGENULL );
  ASSERT ( dopplerPoint->skypos.system == COORDINATESYSTEM_EQUATORIAL,
	   status, FSTATMETRIC_EINPUT, FSTATMETRIC_MSGEINPUT);
  ASSERT ( dopplerPoint->fkdot->length == 2,
	   status, FSTATMETRIC_EINPUT, FSTATMETRIC_MSGEINPUT);

  ASSERT ( multiDetStates, status, FSTATMETRIC_ENULL, FSTATMETRIC_MSGENULL );

  ASSERT ( multidPhi, status, FSTATMETRIC_ENULL, FSTATMETRIC_MSGENULL );
  ASSERT ( *multidPhi == NULL, status, FSTATMETRIC_ENONULL, FSTATMETRIC_MSGENONULL );

  numDet = multiDetStates->length;


  Alpha = dopplerPoint->skypos.longitude;
  Delta = dopplerPoint->skypos.latitude;

  vn[0] = cos(Delta) * cos(Alpha);
  vn[1] = cos(Delta) * sin(Alpha);
  vn[2] = sin(Delta);

  refTimeGPS = multiDetStates->data[0]->data[0].tGPS;	/* use 1st detectors startTime as refTime */
  refTime = GPS2REAL8(refTimeGPS);
  TspanInv = 1.0 / multiDetStates->Tspan;

  /* get tAutumn */
  times.epoch = refTimeGPS;
  TRY ( LALGetEarthTimes ( status->statusPtr, &times ), status );


  if ( (mdPhi = LALCalloc ( 1, sizeof( *mdPhi ))) == NULL ) {
    ABORT (status, FSTATMETRIC_EMEM, FSTATMETRIC_MSGEMEM);
  }
  if ( (mdPhi->data = LALCalloc ( numDet, sizeof( *(mdPhi->data) ) )) == NULL ) {
    ABORT (status, FSTATMETRIC_EMEM, FSTATMETRIC_MSGEMEM);
  }
  mdPhi->length = numDet;

  for ( X=0; X < numDet; X ++ )
    {
      UINT4 numStepsX = multiDetStates->data[X]->length;
      DetectorStateSeries *detStatesX = multiDetStates->data[X];
      PhaseDerivs *dPhi;

      if ( ( dPhi = LALCalloc ( 1, sizeof ( *dPhi ) )) == NULL )
	goto failed;

      dPhi->dFreq = XLALCreateREAL8Vector ( numStepsX );
      dPhi->dAlpha = XLALCreateREAL8Vector ( numStepsX );
      dPhi->dDelta = XLALCreateREAL8Vector ( numStepsX );
      if ( (dPhi->df1dot = XLALCreateREAL8Vector ( numStepsX )) == NULL )
	goto failed;

      mdPhi->data[X] = dPhi;

      for (i=0; i < numStepsX; i++ )
	{
	  REAL8 ti, dT, taui;
	  REAL8 rX[3];	/* radius-vector to use: full, orbital or spin-only */
	  REAL8 rDet[3];	/* vector from earth center to detector */
	  REAL8 fi;		/* instantaneous intrinsic frequency in SSB */
	  PosVel_t posvel;

	  ti = GPS2REAL8 ( detStatesX->data[i].tGPS ) - refTime;

	  /* compute detector-vector relative to earth's center */
	  rDet[0] = detStatesX->data[i].rDetector[0] - detStatesX->data[i].earthState.posNow[0];
	  rDet[1] = detStatesX->data[i].rDetector[1] - detStatesX->data[i].earthState.posNow[1];
	  rDet[2] = detStatesX->data[i].rDetector[2] - detStatesX->data[i].earthState.posNow[2];

	  switch ( phaseType )
	    {
	    case PHASE_FULL:
	      COPY_VECT( rX, detStatesX->data[i].rDetector );
	      break;
	    case PHASE_ORBITAL:
	      COPY_VECT (rX, detStatesX->data[i].earthState.posNow );
	      break;
	    case PHASE_SPIN: /* just given for completeness, probably not very useful */
	      COPY_VECT ( rX, rDet );
	      break;
	    case PHASE_PTOLE: /* use Ptolemaic orbital approximation */
	      getPtolePosVel( &posvel, ti, times.tAutumn );
	      COPY_VECT ( rX, posvel.pos );
	      /* add on the detector-motion due to the Earth's spin */

	      rX[0] += rDet[0];
	      rX[1] += rDet[1];
	      rX[2] += rDet[2];

	      break;
	    default:
	      LALPrintError ("Unknown phase-type specified '%d'\n", phaseType);
	      ABORT ( status, FSTATMETRIC_EINPUT, FSTATMETRIC_MSGEINPUT );
	      break;
	    } /* switch(phaseType) */

	  /* correct for time-delay from SSB to detector */
	  dT = SCALAR(vn, rX );
	  taui = ( ti + dT );

	  fi = dopplerPoint->fkdot->data[0] + taui * dopplerPoint->fkdot->data[1];

	  /* phase-derivatives */
	  dPhi->dFreq->data[i] 	= LAL_TWOPI * offsetUnits->fkdot->data[0] * taui;

	  dPhi->dAlpha->data[i]	= LAL_TWOPI * offsetUnits->skypos.longitude
	    * fi * cos(Delta)*( - rX[0] * sin(Alpha) + rX[1] * cos(Alpha) );

	  dPhi->dDelta->data[i] = LAL_TWOPI * offsetUnits->skypos.latitude
	    * fi * ( - rX[0] * cos(Alpha) * sin(Delta)
		     - rX[1] * sin(Alpha) * sin(Delta)
		     + rX[2] * cos(Delta) );

	  dPhi->df1dot->data[i] = LAL_TWOPI * offsetUnits->fkdot->data[1] * 0.5 * SQ( taui );

	} /* for i < numStepsX */

    } /* for X < numDet */


  (*multidPhi) = mdPhi;

  DETATCHSTATUSPTR (status);
  RETURN(status);

 failed:
  XLALDestroyMultiPhaseDerivs ( mdPhi );
  ABORT ( status, FSTATMETRIC_EMEM, FSTATMETRIC_MSGEMEM );

} /* getMultiPhaseDerivs() */


/** Free all memory */
void
FreeMem ( LALStatus *status, ConfigVariables *cfg )
{

  INITSTATUS (status, "FreeMem", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT ( cfg, status, FSTATMETRIC_ENULL, FSTATMETRIC_MSGENULL );

  TRY ( LALDestroyUserVars ( status->statusPtr ), status );

  LALFree ( cfg->edat->ephemE );
  LALFree ( cfg->edat->ephemS );
  LALFree ( cfg->edat );

  XLALDestroyREAL8Vector ( cfg->dopplerPoint.fkdot );

  gsl_vector_free ( cfg->dopplerOffset );

  XLALDestroyMultiDetectorStateSeries ( cfg->multiDetStates );

  XLALDestroyMultiAMCoeffs ( cfg->multiAMcoe );

  TRY ( LALDestroyMultiNoiseWeights ( status->statusPtr, &(cfg->multiNoiseWeights) ), status );
  XLALDestroyREAL8Vector ( cfg->GLweights );
  XLALDestroyREAL8Vector ( cfg->offsetUnits.fkdot );

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* FreeMem() */


/** Calculate the projected metric onto the subspace of 'c' given by
 * ret_ij = g_ij - ( g_ic * g_jc / g_cc ) , where c is the value of the projected coordinate
 * The output-matrix ret must be allocated
 *
 * return 0 = OK, -1 on error.
 */
int
project_metric( gsl_matrix *ret_ij, gsl_matrix * g_ij, const UINT4 c )
{
  UINT4 i,j;

  if ( !ret_ij )
    return -1;
  if ( !g_ij )
    return -1;
  if ( (ret_ij->size1 != ret_ij->size2) )
    return -1;
  if ( (g_ij->size1 != g_ij->size2) )
    return -1;

  for ( i=0; i < ret_ij->size1; i ++) {
    for ( j=0; j < ret_ij->size2; j ++ ) {
      if ( i==c || j==c ) {
	gsl_matrix_set ( ret_ij, i, j, 0.0 );
      }
      else {
	gsl_matrix_set ( ret_ij, i, j, ( gsl_matrix_get(g_ij, i, j) - (gsl_matrix_get(g_ij, i, c) * gsl_matrix_get(g_ij, j, c) / gsl_matrix_get(g_ij, c, c)) ));
      }
    }
  }
  return 0;
}


/** Calculate the outer product ret_ij of vectors u_i and v_j, given by
 * ret_ij = u_i v_j
 * The output-matrix ret must be allocated and have dimensions len(u) x len(v)
 *
 * return 0 = OK, -1 on error.
 */
int
outer_product (gsl_matrix *ret_ij, const gsl_vector *u_i, const gsl_vector *v_j )
{
  UINT4 i, j;

  if ( !ret_ij || !u_i || !v_j )
    return -1;

  if ( (ret_ij->size1 != u_i->size) || ( ret_ij->size2 != v_j->size) )
    return -1;


  for ( i=0; i < ret_ij->size1; i ++)
    for ( j=0; j < ret_ij->size2; j ++ )
      gsl_matrix_set ( ret_ij, i, j, gsl_vector_get(u_i, i) * gsl_vector_get(v_j, j) );

  return 0;

} /* outer_product() */


/* symmetrize the input-matrix 'mat' (which must be quadratic)
 */
int
symmetrize ( gsl_matrix *mat )
{
  gsl_matrix *tmp;

  if ( !mat )
    return -1;
  if ( mat->size1 != mat->size2 )
    return -1;

  tmp = gsl_matrix_calloc ( mat->size1, mat->size2 );

  gsl_matrix_transpose_memcpy ( tmp, mat );

  gsl_matrix_add (mat, tmp );

  gsl_matrix_scale ( mat, 0.5 );

  gsl_matrix_free ( tmp );

  return 0;

} /* symmetrize() */


/* compute the quadratic form = vec.mat.vec
 */
REAL8
quad_form ( const gsl_matrix *mat, const gsl_vector *vec )
{
  UINT4 i, j;
  REAL8 ret = 0;

  if ( !mat || !vec )
    return 0;
  if ( (mat->size1 != mat->size2) || ( mat->size1 != vec->size ) )
    return 0;

  for (i=0; i < mat->size1; i ++ )
    for (j=0; j < mat->size2; j ++ )
      ret += gsl_vector_get(vec, i) * gsl_matrix_get (mat, i, j) * gsl_vector_get(vec,j);

  return ret;

} /* quad_form() */


/** Get Ptolemaic position and velocity at time tGPS
 * cut-down version of LALDTBaryPtolemaic()
 */

void
getPtolePosVel( PosVel_t *posvel, REAL8 tGPS, REAL8 tAutumnGPS )
{
  REAL8 rev;   /* Earth revolution angle, in radians. */

  /* Some local constants. */
  REAL8 tRev = LAL_AU_SI / LAL_C_SI;
  REAL8 vRev = LAL_TWOPI * tRev / LAL_YRSID_SI;
  REAL8 cosi = cos(LAL_IEARTH);
  REAL8 sini = sin(LAL_IEARTH);

  if ( !posvel )
    return;

  rev = LAL_TWOPI * ( tGPS - tAutumnGPS ) /LAL_YRSID_SI;

  /* Get detector position. */
  posvel->pos[0] = tRev * cos(rev);
  posvel->pos[1] = tRev * sin(rev) * cosi;
  posvel->pos[2]=  tRev * sin(rev) * sini;

  /* Get detector velocity. */
  posvel->vel[0] = -vRev * sin(rev);
  posvel->vel[1] =  vRev * cos(rev) * cosi;
  posvel->vel[2] =  vRev * cos(rev) * sini;

  return;

} /* getPtolePosVel() */

void
XLALDestroyMultiPhaseDerivs ( MultiPhaseDerivs *mdPhi )
{
  UINT4 numDet, X;

  if ( !mdPhi )
    return;
  numDet = mdPhi->length;

  if ( numDet && !mdPhi->data )
    XLAL_ERROR_VOID ( "XLALDestroyMultiPhaseDerivs", XLAL_EINVAL );

  for ( X=0; X < numDet; X ++ )
    {
      PhaseDerivs *dPhiX = mdPhi->data[X];
      if ( !dPhiX )
	continue;
      XLALDestroyREAL8Vector ( dPhiX->dAlpha);
      XLALDestroyREAL8Vector ( dPhiX->dDelta );
      XLALDestroyREAL8Vector ( dPhiX->dFreq );
      XLALDestroyREAL8Vector ( dPhiX->df1dot );

      LALFree ( dPhiX );

    } /* for X < numDet */

  LALFree ( mdPhi->data );
  LALFree ( mdPhi );

  return;

} /* XLALDestroyMultiPhaseDerivs() */

/** Load Ephemeris from ephemeris data-files  */
void
InitEphemeris (LALStatus * status,
	       EphemerisData *edat,	/**< [out] the ephemeris-data */
	       const CHAR *ephemDir,	/**< directory containing ephems */
	       const CHAR *ephemYear,	/**< which years do we need? */
	       LIGOTimeGPS epoch,	/**< epoch of observation */
	       BOOLEAN isLISA		/**< hack this function for LISA ephemeris */
	       )
{
#define FNAME_LENGTH 1024
  CHAR EphemEarth[FNAME_LENGTH];	/* filename of earth-ephemeris data */
  CHAR EphemSun[FNAME_LENGTH];	/* filename of sun-ephemeris data */
  LALLeapSecFormatAndAcc formatAndAcc = {LALLEAPSEC_GPSUTC, LALLEAPSEC_LOOSE};
  INT4 leap;

  INITSTATUS( status, "InitEphemeris", rcsid );
  ATTATCHSTATUSPTR (status);

  ASSERT ( edat, status, FSTATMETRIC_ENULL, FSTATMETRIC_MSGENULL );
  ASSERT ( ephemYear, status, FSTATMETRIC_ENULL, FSTATMETRIC_MSGENULL );

  if ( ephemDir )
    {
      if ( isLISA )
	snprintf(EphemEarth, FNAME_LENGTH, "%s/ephemMLDC.dat", ephemDir);
      else
	snprintf(EphemEarth, FNAME_LENGTH, "%s/earth%s.dat", ephemDir, ephemYear);

      snprintf(EphemSun, FNAME_LENGTH, "%s/sun%s.dat", ephemDir, ephemYear);
    }
  else
    {
      if ( isLISA )
	snprintf(EphemEarth, FNAME_LENGTH, "ephemMLDC.dat");
      else
	snprintf(EphemEarth, FNAME_LENGTH, "earth%s.dat", ephemYear);
      snprintf(EphemSun, FNAME_LENGTH, "sun%s.dat",  ephemYear);
    }

  EphemEarth[FNAME_LENGTH-1]=0;
  EphemSun[FNAME_LENGTH-1]=0;

  /* NOTE: the 'ephiles' are ONLY ever used in LALInitBarycenter, which is
   * why we can use local variables (EphemEarth, EphemSun) to initialize them.
   */
  edat->ephiles.earthEphemeris = EphemEarth;
  edat->ephiles.sunEphemeris = EphemSun;

  TRY (LALInitBarycenter(status->statusPtr, edat), status);

  if ( isLISA )
    {
      edat->leap = -1;	/* dirty hack: signal that ephemeris are in *ECLIPTIC* coords, not EQUATORIAL */
    }
  else
    {
      TRY (LALLeapSecs (status->statusPtr, &leap, &epoch, &formatAndAcc), status);
      edat->leap = (INT2) leap;
    }

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* InitEphemeris() */



/* ================================================================================*/
/* taken from Numerical recipes:
 * function to compute roots and weights for order-N Gauss-Legendre integration
 * modified to use standard C-conventions for arrays (indexed 0... N-1)
 */
#define EPS 3.0e-11

/* Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
 * arrays x[0..n-1] and w[0..n-1] of length n, containing the abscissas and weights of the Gauss-
 * Legendre n-point quadrature formula.
 */
void
gauleg(double x1, double x2, double x[], double w[], int n)
{
     int m,j,i;
     /* High precision is a good idea for this routine. */
     double z1,z,xm,xl,pp,p3,p2,p1;

     /* The roots are symmetric in the interval, so
      * we only have to find half of them. */
     m=(n+1)/2;

     xm=0.5*(x2+x1);
     xl=0.5*(x2-x1);
     /* Loop over the desired roots. */
     for (i=1;i<=m;i++) {
       z=cos( LAL_PI *(i-0.25)/(n+0.5));
       /* Starting with the above approximation to the ith root,
	* we enter the main loop of refinement by Newton's method.
	*/
       do {
	 p1=1.0;
	 p2=0.0;
	 /* Loop up the recurrence relation to get the
	  * Legendre polynomial evaluated at z.
	  */
	 for (j=1;j<=n;j++) {

	   p3=p2;
	   p2=p1;
	   p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
	 }
	 /* p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
	  * by a standard relation involving also p2, the polynomial of one lower order.
	  */
	 pp=n*(z*p1-p2)/(z*z-1.0);
	 z1=z;
	 /* Newton's method. */
	 z=z1-p1/pp;
       } while (fabs(z-z1) > EPS);
       /* Scale the root to the desired interval,  and put in its symmetric counterpart. */
       x[i - 1]=xm-xl*z;	/*RP: changed to C-convention */

       x[n+1-i - 1]=xm+xl*z; 	/*RP: changed to C-convention */

       /* Compute the weight and its symmetric counterpart. */
       w[i-1]=2.0*xl/((1.0-z*z)*pp*pp); /*RP: changed to C-convention */

       w[n+1-i - 1]=w[i - 1];		/*RP: changed to C-convention */
     }

} /* gauleg() */

