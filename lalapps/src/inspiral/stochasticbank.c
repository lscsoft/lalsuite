/*
*  Copyright (C) 2011 Drew Keppel, Duncan Brown
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

/*-----------------------------------------------------------------------
 *
 * File Name: stochasticbank.c
 *
 * Author: Keppel, D.
 *
 *
 *-----------------------------------------------------------------------
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <time.h>
#include <math.h>

#include <lalapps.h>
#include <series.h>
#include <processtable.h>

#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/GenerateInspiral.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Interpolate.h>
#include <lal/TimeSeries.h>
#include <lal/VectorOps.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>
#include "inspiral.h"
#include <LALAppsVCSInfo.h>

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "tmpltbank"

int arg_parse_check( int argc, char *argv[], MetadataTable procparams );
REAL4 ceil_pow_2(REAL4 x);
REAL4 max_mass_ratio(REAL4 mchirp, REAL4 m, REAL4 M);
REAL4FrequencySeries *readPSD(const char *fname, REAL4 fNyq, REAL4 df, UINT4 N, REAL4 fLow, REAL4 fHigh);

typedef struct
tagTemplateWaveformPairs
{
  SnglInspiralTable *sngl_inspiral;
  COMPLEX8FrequencySeries *waveformf;
  struct tagTemplateWaveformPairs *next;
}
TemplateWaveformPairs;

/*
 *
 * variables that control program behaviour
 *
 */

/* debugging */
extern int vrbflg;			/* verbocity of lal function    */

/* template bank generation parameters */
LIGOTimeGPS gpsStartTime = { 0, 0 };	/* input data GPS start time    */
LIGOTimeGPS gpsEndTime	= { 0, 0 };	/* input data GPS end time      */
REAL4   minMass		= -1;		/* minimum component mass       */
REAL4   maxMass		= -1;		/* maximum component mass       */
REAL4   minMTotal	= -1;		/* minimum total mass           */
REAL4   maxMTotal	= -1;		/* maximum total mass           */
REAL4   minMChirp	= -1;		/* minimum chirp mass           */
REAL4   maxMChirp	= -1;		/* maximum chirp mass           */
REAL4   fLo		= -1;		/* low frequency cutoff         */
REAL4   fHi		= -1;		/* high frequency cutoff        */
INT4    facTries	= -1;		/* ratio of numTries/numTmplts  */
REAL4   minMatch	= -1;		/* minimal match of tmplt bank  */
enum { unset, urandom, user } randSeedType = unset;    /* sim seed type */
INT4    randomSeed	= 0;		/* value of sim rand seed       */
LIGOTimeGPS tc		= LIGOTIMEGPSZERO;

/* input/output parameters */
CHAR   *userTag		= NULL;
CHAR   *psdFileName	= NULL;
INT4   outCompress	= 0;

REAL4 ceil_pow_2(REAL4 x)
{
  INT4 exponent = 0;
  REAL4 mantissa;

  mantissa = frexp(x, &exponent);

  if (mantissa < 0)
    return NAN;
  else if (mantissa == 0.5)
    return ldexp((REAL4) 1., exponent - 1);
  else
    return ldexp((REAL4) 1., exponent);
}


REAL4 max_mass_ratio(REAL4 mchirp, REAL4 m, REAL4 M)
{
  REAL4 qmax1,qmax2,tmp1,tmp2;

  tmp1 = pow(mchirp / m, 5.);
  tmp2 = tmp1 * (-9. + pow(3., .5) * pow(27. + 4.*tmp1, .5) ) / 18.;
  tmp2 = pow(tmp2, 1./3.);
  qmax1 = tmp2 + tmp1 / 3. / tmp2;

  tmp1 = -M;
  tmp2 = pow(mchirp, 5./3.) * pow(M, 1./3.);
  qmax2 = -tmp1 + sqrt(tmp1*tmp1 - 4.*tmp2);
  qmax2 /= 2. * m;

  return qmax1 < qmax2 ? qmax1 : qmax2;
}


REAL4FrequencySeries *readPSD(const char *fname, REAL4 fNyq, REAL4 df, UINT4 N, REAL4 fLow, REAL4 fHigh)
{
  FILE *fp;
  UINT4 i, j;
  REAL4 f0, f1;
  REAL8 junk;
  char line[1000];
  REAL8FrequencySeries *psd1;
  REAL4FrequencySeries *psd;

  fp = fopen(fname, "r");
  i = 0;
  while (fgets(line, sizeof line, fp))
  {
    if (*line == '#')
      continue;
    else if (i == 0)
      sscanf(line, "%e\t%le", &f0, &junk);
    else if (i == 1)
      sscanf(line, "%e\t%le", &f1, &junk);
    i++;
  }
  fclose(fp);

  psd1 = XLALCreateREAL8FrequencySeries("", &tc, f0, f1 - f0, &lalSecondUnit, i);
  fp = fopen(fname, "r");
  i = 0;
  while (fgets(line, sizeof line, fp))
  {
    if (*line == '#')
      continue;
    sscanf(line, "%le\t%le", &junk, &(psd1->data->data[i]));
    i++;
  }
  fclose(fp);


  REAL8 Fs[psd1->data->length];
  for (i = 0; i < psd1->data->length; i++)
  {
    Fs[i] = (REAL8) (i * (REAL4) psd1->deltaF + (REAL4) psd1->f0);
  }

  psd = XLALCreateREAL4FrequencySeries("", &tc, -fNyq, df, &lalSecondUnit, N);

  for (i = 0; i < N / 2 + 1; i++)
  {
    REAL8 tmpS, dS;
    REAL8 tmpf = (REAL8) fNyq - df*i;

    if ( tmpf < fLow || tmpf > fHigh )
    {
      psd->data->data[i] = 0.;
      if (i > 0)
        psd->data->data[N - i] = 0.;
      continue;
    }

    for (j=2; j < 10; j++)
    {
      dS = XLALREAL8PolynomialInterpolation(&tmpS, tmpf, psd1->data->data, Fs, 10);
      if ( (tmpS == 0. && dS < 1e-5) || (fabs(dS / tmpS) < 1e-5) )
        break;
    }

    psd->data->data[i] = tmpS;
    if (i > 0)
      psd->data->data[N - i] = tmpS;
  }
  XLALDestroyREAL8FrequencySeries(psd1);

  return psd;
}


int main ( int argc, char *argv[] )
{
  /* lal function variables */
  LALStatus             status = blank_status;

  /* templates */
  RandomParams         *randParams = NULL;
  InspiralTemplate      newTmplt;
  SnglInspiralTable    *thisTmplt  = NULL;

  /* output data */
  MetadataTable         templateBank;
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         searchsummvars;
  SearchSummvarsTable  *this_search_summvar = NULL;
  ProcessParamsTable   *this_proc_param = NULL;
  LIGOLwXMLStream       results;

  /* waveform storage */

  REAL4TimeSeries *realWaveform = NULL;		/* storing real waveform */
  REAL4TimeSeries *imagWaveform = NULL;		/* storing imag waveform */
  REAL4TimeSeries *overlap = NULL;		/* storing abs(complex overlap) */
  COMPLEX8TimeSeries *waveform = NULL;		/* storing complex waveform */
  COMPLEX8TimeSeries *overlapc = NULL;		/* storing complex overlap */

  COMPLEX8FrequencySeries *waveformf1 = NULL;	/* storing FFT of complex waveform 1 */
  COMPLEX8FrequencySeries *waveformf2 = NULL;	/* storing FFT of complex waveform 2 */
  COMPLEX8FrequencySeries *overlapf = NULL;	/* storing FFT of complex overlap */

  REAL4FrequencySeries *psd = NULL;		/* storing inverse psd */

  TemplateWaveformPairs *headWvTmpltPr = NULL;
  TemplateWaveformPairs *thisWvTmpltPr = NULL;
  TemplateWaveformPairs *tmpWvTmpltPr = NULL;

  /* counters and other variables */
  UINT4 i, N, trynum, lastTmplt, numTmplts;
  INT4 fs;
  CHAR fname[256];
  REAL4 dt, df, norm, match;
  COMPLEX8FFTPlan *fwdp = NULL;
  COMPLEX8FFTPlan *revp = NULL;

  /*
   *
   * initialization
   *
   */


  /* set up inital debugging values */
  XLALSetErrorHandler(XLALAbortErrorHandler);
  lal_errhandler = LAL_ERR_EXIT;

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *)
    calloc( 1, sizeof(ProcessTable) );
  XLALGPSTimeNow(&(proctable.processTable->start_time));
  XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME, lalAppsVCSIdentId,
      lalAppsVCSIdentStatus, lalAppsVCSIdentDate, 0);
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *)
    calloc( 1, sizeof(ProcessParamsTable) );

  /* call the argument parse and check function */
  arg_parse_check( argc, argv, procparams );

  /* can use LALMalloc() / LALCalloc() from here */


  /*
   *
   * create the radom number seed
   *
   */


  /* store the seed in the search summvars table */
  this_search_summvar = searchsummvars.searchSummvarsTable =
    (SearchSummvarsTable *) LALCalloc( 1, sizeof(SearchSummvarsTable) );
  snprintf( this_search_summvar->name,
      LIGOMETA_NAME_MAX, "template bank simulation seed" );

  if ( randSeedType == urandom )
  {
    FILE   *fpRand = NULL;
    INT4    randByte;

    if ( vrbflg )
      fprintf( stdout, "obtaining random seed from /dev/urandom: " );

    randomSeed = 0;
    fpRand = fopen( "/dev/urandom", "r" );
    if ( fpRand )
    {
      for ( randByte = 0; randByte < 4 ; ++randByte )
      {
        INT4 tmpSeed = (INT4) fgetc( fpRand );
        randomSeed += tmpSeed << ( randByte * 8 );
      }
      fclose( fpRand );
    }
    else
    {
      perror( "error obtaining random seed from /dev/urandom" );
      exit( 1 );
    }
  }
  else if ( randSeedType == user )
  {
    if ( vrbflg )
      fprintf( stdout, "using user specified random seed: " );
  }

  this_search_summvar->value = randomSeed;
  snprintf( this_search_summvar->string, LIGOMETA_STRING_MAX,
      "%d", randomSeed );
  if ( vrbflg ) fprintf( stdout, "%d\n", randomSeed );

  /* create the tmplt bank random parameter structure */
  LAL_CALL( LALCreateRandomParams( &status, &randParams, randomSeed ),
      &status );

  // compute length of longest waveform
  memset( &newTmplt, 0, sizeof(InspiralTemplate) );
  newTmplt.massChoice = m1Andm2;
  newTmplt.order = LAL_PNORDER_TWO;
  newTmplt.fLower = fLo;
  newTmplt.fCutoff = fHi;
  newTmplt.mass1 = minMChirp / pow(.25, .6) / 2;
  newTmplt.mass2 = newTmplt.mass1;

  LAL_CALL( LALInspiralParameterCalc( &status, &newTmplt ), &status );

  fs = 2 * ceil_pow_2(fHi * 1.1);
  dt = 1. / fs;
  N = fs * 2 * ceil_pow_2(newTmplt.t0 * 1.1 + 32.);
  df = (REAL4) fs / N;

  memset( &newTmplt, 0, sizeof(InspiralTemplate) );
  newTmplt.massChoice = m1Andm2;
  newTmplt.approximant = TaylorT1;
  newTmplt.order = LAL_PNORDER_THREE_POINT_FIVE;
  newTmplt.ampOrder = LAL_PNORDER_NEWTONIAN;
  newTmplt.fLower = fLo;
  newTmplt.fCutoff = fHi;
  newTmplt.tSampling = fs;
  newTmplt.signalAmplitude = 1.;

  fwdp = XLALCreateForwardCOMPLEX8FFTPlan(N, 1);
  revp = XLALCreateReverseCOMPLEX8FFTPlan(N, 1);

  psd = readPSD(psdFileName, fs/2, df, N, fLo, fHi);

  realWaveform = XLALCreateREAL4TimeSeries("", &tc, 0.0, dt, &lalDimensionlessUnit, N);
  imagWaveform = XLALCreateREAL4TimeSeries("", &tc, 0.0, dt, &lalDimensionlessUnit, N);
  waveform = XLALCreateCOMPLEX8TimeSeries("", &tc, 0.0, dt, &lalDimensionlessUnit, N);
  overlapc = XLALCreateCOMPLEX8TimeSeries("", &tc, 0.0, dt, &lalDimensionlessUnit, N);
  overlap = XLALCreateREAL4TimeSeries("", &tc, 0.0, dt, &lalDimensionlessUnit, N);

  waveformf1 = XLALCreateCOMPLEX8FrequencySeries("", &tc, 0.0, df, &lalSecondUnit, N);
  overlapf = XLALCreateCOMPLEX8FrequencySeries("", &tc, 0.0, df, &lalSecondUnit, N);

  /*
   *
   * create a stochastic template bank
   *
   */


  /* make sure the pointer to the first template is null */
  
  templateBank.snglInspiralTable = NULL;

  trynum = 0;
  numTmplts = 0;
  lastTmplt = 0;
  while ( trynum <= facTries * numTmplts )
  {
    fprintf(stderr, "\r%i %i %i", trynum, numTmplts, lastTmplt);

    REAL4 randfloat,mchirp,q,m1,m2,M,eta;
    REAL4 qmax;

    trynum++;
    lastTmplt++;

    /* generate random parameters for the injection */
    randfloat = XLALUniformDeviate( randParams );
    if ( XLAL_IS_REAL4_FAIL_NAN( randfloat ) )
    {
      exit( 1 );
    }

    mchirp = (maxMChirp - minMChirp) * randfloat;
    mchirp += minMChirp;

    qmax = max_mass_ratio(mchirp, minMass, maxMTotal);

    m2 = 0;
    m1 = 0;
    while (m2 < minMass)
    {
      /* generate random parameters for the injection */
      randfloat = XLALUniformDeviate( randParams );
      if ( XLAL_IS_REAL4_FAIL_NAN( randfloat ) )
      {
        exit( 1 );
      }

      q = (qmax - 1.) * randfloat;
      q += 1.;

      eta = q;
      eta /= 1.+q;
      eta /= 1.+q;

      M = mchirp / pow(eta, 0.6);
      m2 = M / (q + 1.);
      m1 = M - m2;
    }

    newTmplt.mass1 = m1;
    newTmplt.mass2 = m2;


    LAL_CALL( LALInspiralParameterCalc( &status, &newTmplt ), &status );


    /* --- now we can call the injection function --- */
    memset(realWaveform->data->data, 0., realWaveform->data->length*sizeof(REAL4));
    memset(imagWaveform->data->data, 0., imagWaveform->data->length*sizeof(REAL4));
    memset(waveform->data->data, 0., waveform->data->length*sizeof(COMPLEX8));

    LAL_CALL(LALInspiralWaveTemplates(&status, realWaveform->data, imagWaveform->data, &newTmplt), &status);

    for (i=0; i< realWaveform->data->length; i++)
    {
      waveform->data->data[i] = crectf( realWaveform->data->data[i], imagWaveform->data->data[i] );
    }

    // FFT complex waveform

    XLALCOMPLEX8TimeFreqFFT(waveformf1, waveform, fwdp);

    // whiten waveform

    XLALWhitenCOMPLEX8FrequencySeries(waveformf1, psd);

    // normalize waveform

    XLALCCVectorMultiplyConjugate(overlapf->data, waveformf1->data, waveformf1->data);
    XLALCOMPLEX8FreqTimeFFT(overlapc, overlapf, revp);
    XLALCOMPLEX8VectorAbs(overlap->data, overlapc->data);
    norm = 0.;
    for (i = 0; i < overlap->data->length; i++)
    {
      if (overlap->data->data[i] > norm)
        norm = overlap->data->data[i];
    }
    norm = sqrt(norm);

    for (i = 0; i < waveformf1->data->length; i++)
    {
      waveformf1->data->data[i] /= norm;
    }

    // compute overlap with previous tmplts

    tmpWvTmpltPr = headWvTmpltPr;
    match = 0.;
    while (tmpWvTmpltPr)
    {
      waveformf2 = tmpWvTmpltPr->waveformf;
      XLALCCVectorMultiplyConjugate(overlapf->data, waveformf1->data, waveformf2->data);
      XLALCOMPLEX8FreqTimeFFT(overlapc, overlapf, revp);
      XLALCOMPLEX8VectorAbs(overlap->data, overlapc->data);
      match = 0.;
      for (i = 0; i < overlap->data->length; i++)
      {
        if (overlap->data->data[i] > match)
          match = overlap->data->data[i];
      }

      if (match > minMatch)
        break;
      tmpWvTmpltPr = tmpWvTmpltPr->next;
    }

    if (match > minMatch)
      continue;

    // save waveform in tmplt list

    if ( ! templateBank.snglInspiralTable )
    {
      thisTmplt = templateBank.snglInspiralTable =
        (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
    }
    else
    {
      thisTmplt = thisTmplt->next =
        (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
    }

    thisTmplt->mass1 = newTmplt.mass1;
    thisTmplt->mass2 = newTmplt.mass2;
    thisTmplt->mchirp = newTmplt.chirpMass;
    thisTmplt->eta = newTmplt.eta;
    thisTmplt->tau0 = newTmplt.t0;
    thisTmplt->tau2 = newTmplt.t2;
    thisTmplt->tau3 = newTmplt.t3;
    thisTmplt->tau4 = newTmplt.t4;
    thisTmplt->tau5 = newTmplt.t5;
    thisTmplt->ttotal = newTmplt.tC;
    thisTmplt->psi0 = newTmplt.psi0;
    thisTmplt->psi3 = newTmplt.psi3;
    thisTmplt->f_final = newTmplt.fFinal;
    thisTmplt->eta = newTmplt.eta;
    thisTmplt->beta = newTmplt.beta;
    snprintf( thisTmplt->ifo, LIGOMETA_IFO_MAX, "P1" );
    snprintf( thisTmplt->search, LIGOMETA_SEARCH_MAX, "stochasticbank" );
    snprintf( thisTmplt->channel, LIGOMETA_CHANNEL_MAX, "SIM-BANK" );

    if (thisWvTmpltPr == NULL)
      thisWvTmpltPr = (TemplateWaveformPairs *) LALCalloc(1, sizeof(TemplateWaveformPairs));
    else
    {
      thisWvTmpltPr->next = (TemplateWaveformPairs *) LALCalloc(1, sizeof(TemplateWaveformPairs));
      thisWvTmpltPr = thisWvTmpltPr->next;
    }
    thisWvTmpltPr->sngl_inspiral = thisTmplt;
    thisWvTmpltPr->waveformf = waveformf1;
    if (headWvTmpltPr == NULL)
      headWvTmpltPr = thisWvTmpltPr;

    waveformf1 = XLALCreateCOMPLEX8FrequencySeries("", &tc, 0.0, df, &lalSecondUnit, N);

    numTmplts++;
    lastTmplt = 0;
  }
  fprintf(stderr, "\r%i %i %i\n", trynum, numTmplts, lastTmplt);

  /*
   *
   * write the output data
   *
   */

  /* open the output xml file */
  memset( &results, 0, sizeof(LIGOLwXMLStream) );
  if ( userTag && !outCompress )
  {
    snprintf( fname, sizeof(fname), "P1-TMPLTBANK_%s-%d-%d.xml",
        userTag, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else if ( userTag && outCompress )
  {
    snprintf( fname, sizeof(fname), "P1-TMPLTBANK_%s-%d-%d.xml.gz",
        userTag, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else if ( !userTag && outCompress )
  {
    snprintf( fname, sizeof(fname), "P1-TMPLTBANK-%d-%d.xml.gz",
        gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else
  {
    snprintf( fname, sizeof(fname), "P1-TMPLTBANK-%d-%d.xml",
        gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &results, fname ), &status );

  /* write the process table */
  snprintf( proctable.processTable->ifos, LIGOMETA_IFO_MAX, "P1" );
  XLALGPSTimeNow(&(proctable.processTable->end_time));
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_table ),
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, proctable,
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  free( proctable.processTable );

  /* erase the first empty process params entry */
  {
    ProcessParamsTable *emptyPPtable = procparams.processParamsTable;
    procparams.processParamsTable = procparams.processParamsTable->next;
    free( emptyPPtable );
  }

  /* write the process params table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_params_table ),
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, procparams,
        process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  while( procparams.processParamsTable )
  {
    this_proc_param = procparams.processParamsTable;
    procparams.processParamsTable = this_proc_param->next;
    free( this_proc_param );
  }

  /* write the search summvars table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results,
        search_summvars_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, searchsummvars,
        search_summvars_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  while( searchsummvars.searchSummvarsTable )
  {
    this_search_summvar = searchsummvars.searchSummvarsTable;
    searchsummvars.searchSummvarsTable = this_search_summvar->next;
    LALFree( this_search_summvar );
  }

  /* write the template bank to the file */
  if ( templateBank.snglInspiralTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, sngl_inspiral_table ),
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, templateBank,
          sngl_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  }

  while ( templateBank.snglInspiralTable )
  {
    thisTmplt = templateBank.snglInspiralTable;
    templateBank.snglInspiralTable = templateBank.snglInspiralTable->next;
    LALFree( thisTmplt );
  }

  /* close the output xml file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &results ), &status );

  /* free the rest of the memory, check for memory leaks and exit */

  LAL_CALL( LALDestroyRandomParams( &status, &randParams ), &status );

  while (headWvTmpltPr)
  {
    thisWvTmpltPr = headWvTmpltPr;
    headWvTmpltPr = headWvTmpltPr->next;
    XLALDestroyCOMPLEX8FrequencySeries(thisWvTmpltPr->waveformf);
    LALFree(thisWvTmpltPr);
  }

  XLALDestroyCOMPLEX8FFTPlan(fwdp);
  XLALDestroyCOMPLEX8FFTPlan(revp);
  XLALDestroyREAL4TimeSeries(realWaveform);
  XLALDestroyREAL4TimeSeries(imagWaveform);
  XLALDestroyREAL4TimeSeries(overlap);
  XLALDestroyREAL4FrequencySeries(psd);
  XLALDestroyCOMPLEX8TimeSeries(waveform);
  XLALDestroyCOMPLEX8TimeSeries(overlapc);
  XLALDestroyCOMPLEX8FrequencySeries(waveformf1);
  XLALDestroyCOMPLEX8FrequencySeries(overlapf);
  LALCheckMemoryLeaks();
  exit( 0 );
}


/* ------------------------------------------------------------------------- */

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  calloc( 1, sizeof(ProcessParamsTable) ); \
  snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
      PROGRAM_NAME ); \
      snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
          long_options[option_index].name ); \
          snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
          snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

#define USAGE \
"  --help                       display this message\n"\
"  --version                    print version information and exit\n"\
"  --user-tag STRING            set the process_params usertag to STRING\n"\
"\n"\
"  --gps-start-time SEC         GPS second of data start time\n"\
"  --gps-end-time SEC           GPS second of data end time\n"\
"  --low-frequency-cutoff F     Compute tau parameters from F Hz\n"\
"  --minimum-mass MASS          set minimum component mass of bank to MASS\n"\
"  --maximum-mass MASS          set maximum component mass of bank to MASS\n"\
"  --tries-factor N             test a factor N more points than templates retained\n"\
"  --random-seed SEED           set random number seed for injections to SEED\n"\
  "                                 (urandom|integer)\n"\
"  --write-compress             write a compressed xml file\n"\
  "\n"

int arg_parse_check( int argc, char *argv[], MetadataTable procparams )
{
  /* getopt arguments */
  struct option long_options[] =
  {
    /* these options set a flag */
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"write-compress",          no_argument,       &outCompress,      1 },
    {"help",                    no_argument,       0,                'h'},
    {"gps-start-time",          required_argument, 0,                's'},
    {"gps-end-time",            required_argument, 0,                'e'},
    {"low-frequency-cutoff",    required_argument, 0,                'f'},
    {"high-frequency-cutoff",   required_argument, 0,                'F'},
    {"minimum-mass",            required_argument, 0,                'a'},
    {"maximum-mass",            required_argument, 0,                'A'},
    {"minimum-mtotal",          required_argument, 0,                'b'},
    {"maximum-mtotal",          required_argument, 0,                'B'},
    {"minimum-mchirp",          required_argument, 0,                'c'},
    {"maximum-mchirp",          required_argument, 0,                'C'},
    {"tries-factor",            required_argument, 0,                'N'},
    {"minimal-match",           required_argument, 0,                'm'},
    {"psd-file",                required_argument, 0,                'P'},
    {"random-seed",             required_argument, 0,                'S'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {0, 0, 0, 0}
  };
  int c;
  ProcessParamsTable *this_proc_param = procparams.processParamsTable;


  /*
   *
   * parse command line arguments
   *
   */


  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, "ha:b:c:e:m:s:A:B:C:N:P:S:Z:",
        long_options, &option_index );

    /* detect the end of the options */
    if ( c == - 1 )
    {
      break;
    }

    switch ( c )
    {
      case 0:
        /* if this option set a flag, do nothing else now */
        if ( long_options[option_index].flag != 0 )
        {
          break;
        }
        else
        {
          fprintf( stderr, "error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;

      case 'h':
        fprintf( stdout, USAGE );
        exit( 0 );
        break;

      case 'Z':
        /* create storage for the usertag */
        optarg_len = strlen( optarg ) + 1;
        userTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( userTag, optarg, optarg_len );

        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s",
            PROGRAM_NAME );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--userTag" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            optarg );
        break;

      case 'P':
        /* create storage for the psd file name */
        optarg_len = strlen( optarg ) + 1;
        psdFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( psdFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", psdFileName );
        break;

      case 'a':
        minMass = (REAL4) atof( optarg );
        if ( minMass <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "miniumum component mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, minMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", minMass );
        break;

      case 'A':
        maxMass = (REAL4) atof( optarg );
        if ( maxMass <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "maxiumum component mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, maxMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", maxMass );
        break;

      case 'b':
        minMTotal = (REAL4) atof( optarg );
        if ( minMTotal <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "miniumum total mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, minMTotal );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", minMTotal );
        break;

      case 'B':
        maxMTotal = (REAL4) atof( optarg );
        if ( maxMTotal <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "maxiumum total mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, maxMTotal );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", maxMTotal );
        break;

      case 'c':
        minMChirp = (REAL4) atof( optarg );
        if ( minMChirp <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "miniumum chirp mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, minMChirp );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", minMChirp );
        break;

      case 'C':
        maxMChirp = (REAL4) atof( optarg );
        if ( maxMChirp <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "maxiumum chirp mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, maxMChirp );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", maxMChirp );
        break;

      case 'f':
        fLo = (REAL4) atof( optarg );
        if ( fLo <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "lower frequency cutoff must be > 0: "
              "(%f Hz specified)\n",
              long_options[option_index].name, fLo );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", fLo );
        break;

      case 'F':
        fHi = (REAL4) atof( optarg );
        if ( fHi <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "upper frequency cutoff must be > 0: "
              "(%f Hz specified)\n",
              long_options[option_index].name, fHi );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", fHi );
        break;

      case 'J':
        if ( ! strcmp( "urandom", optarg ) )
        {
          randSeedType = urandom;
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        else
        {
          randSeedType = user;
          randomSeed = (INT4) atoi( optarg );
          ADD_PROCESS_PARAM( "int", "%d", randomSeed );
        }
        break;

      case 'N':
        facTries = (INT4) atoi( optarg );
        if ( facTries < 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "ratio of test points to templates"
              "must be greater than 1: (%d specified)\n",
              long_options[option_index].name, facTries );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", facTries );
        break;

      case 'm':
        minMatch = (REAL4) atof( optarg );
        if ( minMatch < 0 || minMatch >= 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "minimal match of the template bank"
              "must be in [0,1): (%f specified)\n",
              long_options[option_index].name, minMatch );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%f", minMatch );
        break;

      case 's':
        {
          long int gstartt = atol( optarg );
          if ( gstartt < 441417609 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS start time is prior to "
                "Jan 01, 1994  00:00:00 UTC:\n"
                "(%ld specified)\n",
                long_options[option_index].name, gstartt );
            exit( 1 );
          }
          if ( gstartt > 999999999 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS start time is after "
                "Sep 14, 2011  01:46:26 UTC:\n"
                "(%ld specified)\n",
                long_options[option_index].name, gstartt );
            exit( 1 );
          }
          gpsStartTime.gpsSeconds = (INT4) gstartt;
          gpsStartTime.gpsNanoSeconds = 0;
          ADD_PROCESS_PARAM( "int", "%ld", gstartt );
        }
        break;

      case 'e':
        {
          long int gendt = atol( optarg );
          if ( gendt > 999999999 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS end time is after "
                "Sep 14, 2011  01:46:26 UTC:\n"
                "(%ld specified)\n",
                long_options[option_index].name, gendt );
            exit( 1 );
          }
          else if ( gendt < 441417609 )
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "GPS end time is prior to "
                "Jan 01, 1994  00:00:00 UTC:\n"
                "(%ld specified)\n",
                long_options[option_index].name, gendt );
            exit( 1 );
          }
          gpsEndTime.gpsSeconds = (INT4) gendt;
          gpsEndTime.gpsNanoSeconds = 0;
          ADD_PROCESS_PARAM( "int", "%ld", gendt );
        }
        break;

      case '?':
        fprintf( stderr, USAGE );
        exit( 1 );
        break;

      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        fprintf( stderr, USAGE );
        exit( 1 );
    }
  }

  if ( optind < argc )
  {
    fprintf( stderr, "extraneous command line arguments:\n" );
    while ( optind < argc )
    {
      fprintf ( stderr, "%s\n", argv[optind++] );
    }
    exit( 1 );
  }

  /*
   *
   * check validity of arguments
   *
   */


  if ( minMass < 0 )
  {
    fprintf( stderr, "--minimum-mass must be specified\n" );
    exit( 1 );
  }
  if ( maxMass < 0 )
  {
    fprintf( stderr, "--maximum-mass must be specified\n" );
    exit( 1 );
  }

  if ( minMTotal < 0 )
  {
    fprintf( stderr, "--minimum-mtotal must be specified\n" );
    exit( 1 );
  }
  if ( maxMTotal < 0 )
  {
    fprintf( stderr, "--maximum-mtotal must be specified\n" );
    exit( 1 );
  }

  if ( minMChirp < 0 )
  {
    fprintf( stderr, "--minimum-mchirp must be specified\n" );
    exit( 1 );
  }
  if ( maxMChirp < 0 )
  {
    fprintf( stderr, "--maximum-mchirp must be specified\n" );
    exit( 1 );
  }

  if ( fLo < 0 )
  {
    fprintf( stderr, "--low-frequency-cutoff must be specified\n" );
    exit( 1 );
  }
  if ( fHi < 0 )
  {
    fprintf( stderr, "--high-frequency-cutoff must be specified\n" );
    exit( 1 );
  }

  /* check that the bank parameters have been specified */
  if ( facTries < 0 )
  {
    fprintf( stderr, "--tries-factor must be specified\n" );
    exit( 1 );
  }

  return 0;
}

#undef ADD_PROCESS_PARAM
