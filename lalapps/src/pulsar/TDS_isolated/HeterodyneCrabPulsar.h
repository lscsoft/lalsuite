/* Matt Pitkin 10/03/04 - TNHeterodyne struct changed */

/* Matt Pitkin 26/03/04 - function arguments have been changed to conform to the
   LAL Spec */
/*
$Id$
*/

#ifndef _HETERODYNECRABPULSAR_H
#define _HETERODYNECRABPULSAR_H

#include <lal/LALStdlib.h>

#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/SkyCoordinates.h>
#include <lal/AVFactories.h>

#ifdef __cplusplus
extern "C" {
#endif

/* DEFINE RCS ID STRING */
NRCSID (HETERODYNECRABPULSARH, "$Id$"); 

/* <lalErrTable file="GetCrabParamsErrorTable"> */
#define HETERODYNECRABPULSARH_ENULLINPUT 1
#define HETERODYNECRABPULSARH_ENULLOUTPUT 2
#define HETERODYNECRABPULSARH_ENULLVALUES 3
#define HETERODYNECRABPULSARH_EEPHEMERISFILENAME 4
#define HETERODYNECRABPULSARH_ENUMEPHEMERISDATA 5
#define HETERODYNECRABPULSARH_ELENGTHINPUTOUTPUT 6
#define HETERODYNECRABPULSARH_EINVALIDF0 7

#define HETERODYNECRABPULSARH_MSGENULLINPUT "Input was Null"
#define HETERODYNECRABPULSARH_MSGENULLOUTPUT "Output was Null"
#define HETERODYNECRABPULSARH_MSGENULLVALUES "Values was Null"
#define HETERODYNECRABPULSARH_MSGEEPHEMERISFILENAME "No ephemeris filename given"
#define HETERODYNECRABPULSARH_MSGENUMEPHEMERISDATA "Number of ephemeris data points must be greater than zero"
#define HETERODYNECRABPULSARH_MSGELENGTHINPUTOUTPUT "Input and output vectors are not the right lengths"
#define HETERODYNECRABPULSARH_MSGEINVALIDF0 "Input value of F0 is invalid, must be greater than 0"
/* </lalErrTable> */

/* define new structures and type for crab ephemeris reading - to be put in header file later */
typedef struct
tagGetCrabEphemerisInput
{
	CHAR *filename; /* file name of the Crab ephemeris (crab_ephemeris.txt) */
} GetCrabEphemerisInput;

typedef struct
tagCrabSpindownParamsInput
{
  REAL8Vector *tArr;	/* vector of pulse arrival times from Crab ephemeris (GPS seconds) */
  REAL8Vector *f0;	/* vector of pulsar frequencies at each t_arr */
  REAL8Vector *f1;	/* vector of first frequecy derivs at each t_arr */
	UINT4 numOfData;  /* num of data in ephemeris file */
} CrabSpindownParamsInput;

typedef struct
tagCrabSpindownParamsOutput
{
  REAL8Vector *tArr;	/* vector of pulse arrival times from Crab ephemeris (GPS seconds) */
  REAL8Vector *f0;	/* vector of pulsar frequencies at each t_arr */
  REAL8Vector *f1;	/* vector of first frequecy derivs at each t_arr */
  REAL8Vector *f2;	/* vectors of second, third and fourth frequency */
  REAL8Vector *f3;	/* derivs of the Crab pulsar calculated at each  */
  REAL8Vector *f4;	/* t_arr, length will be one less than t_arr in input    */
	UINT4 numOfData;  /* num of data in ephemeris file */
} CrabSpindownParamsOutput;

/*typedef struct
tagParamsForHeterodyne
{
  REAL8 f0;
  REAL8 f1;
  REAL8 f2;
  REAL8 f3;
  REAL8 f4;
  REAL8 epoch;
  UINT4 startPos;
  EphemerisData *edat;
  LALDetector detector;
} ParamsForHeterodyne;*/

/* ParamsForHeterodyne struct changed */
typedef struct
tagParamsForHeterodyne
{
  REAL8 f0;
  REAL8 f1;
  REAL8 f2;
  REAL8 f3;
  REAL8 f4;
  REAL8 epoch;
} ParamsForHeterodyne;

/*typedef struct
tagTNHeterodyneInput
{
  REAL8			f0;
  REAL8			f1;
  REAL8			f2;
  COMPLEX16TimeSeries	Vh;
  REAL8			t0;
  REAL8			phase;
  SkyPosition           source;  
} TNHeterodyneInput;*/

/*TNHeterodyne structure that doesn't need to heterodyne at the SSB */
typedef struct
tagTNHeterodyneInput
{
  REAL8			f0;
  REAL8			f1;
  REAL8			f2;
  COMPLEX16 Vh;
	REAL8			t0;
	LIGOTimeGPS epoch;  
} TNHeterodyneInput; 

typedef struct
tagTNHeterodyneOutput
{
  COMPLEX16	Vh;
	REAL8 Dphase;
} TNHeterodyneOutput;

/* function definitions */

void 
LALGetCrabEphemeris	( LALStatus			*status,
			  CrabSpindownParamsInput	*output,
				GetCrabEphemerisInput *input );
			  
void
LALComputeFreqDerivatives	( LALStatus			*status,
				  CrabSpindownParamsOutput	*ouput,
					CrabSpindownParamsInput	*input );			  

void
LALSetSpindownParams	( LALStatus			*status,
			  ParamsForHeterodyne		*output,
				CrabSpindownParamsOutput	*input,
				LIGOTimeGPS			epoch );

void
LALTimingNoiseHeterodyne	( LALStatus		*status,
				  TNHeterodyneOutput	*output,
					TNHeterodyneInput	*input,
				  ParamsForHeterodyne	*params );

#ifdef __cplusplus
}
#endif
				  
#endif /* _HETERODYNECRABPULSAR_H */
