/*
*  Copyright (C) 2007 John Whelan
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

/* ---------- see CoarseGrainFrequencySeries.h for doxygen documentation ---------- */

#include <complex.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/CoarseGrainFrequencySeries.h>
#include <math.h>


void
LALSCoarseGrainFrequencySeries(LALStatus                      *status,
                               REAL4FrequencySeries           *output,
                               const REAL4FrequencySeries     *input,
                               const FrequencySamplingParams  *params)

{
  UINT4         lengthCoarse, lengthFine;
  REAL8         f0Coarse, f0Fine;
  REAL8         fMinCoarse, fMinFine;
  REAL8         deltaFCoarse, deltaFFine;
  REAL8         offset, resRatio;

  UINT4         k, l;
  UINT4         lMin, lMax;
  REAL8         lamMin, lamMax;
  REAL4         value;

  /* initialize status structure */
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  /* checks for null pointers: */

  /*    output series */
  ASSERT(output != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of output series */
  ASSERT(output->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of data member of output series */
  ASSERT(output->data->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    input series */
  ASSERT(input != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of input series */
  ASSERT(input->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of data member of input series */
  ASSERT(input->data->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    parameter structure */
  ASSERT(params != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /* checks for duplicate pointers: */
  ASSERT(output != input, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  ASSERT(output->data != input->data, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  ASSERT(output->data->data != input->data->data, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  /* extract coarse-grained parameters  */
  lengthCoarse = params->length;
  f0Coarse     = params->f0;
  deltaFCoarse = params->deltaF;
  fMinCoarse = ( f0Coarse == 0.0 ? 0.0 : f0Coarse - deltaFCoarse / 2.0 );

  /* extract fine-grained parameters */
  lengthFine   = input->data->length;
  f0Fine       = input->f0;
  deltaFFine   = input->deltaF;
  fMinFine = ( f0Fine == 0.0 ? 0.0 : f0Fine - deltaFFine / 2.0 );

  /* check for legality of values */

  /*    length must be positive */

  ASSERT(lengthCoarse != 0, status,
         COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

  ASSERT(lengthFine != 0, status,
         COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

  /*    start frequency must not be negative */

  if (fMinCoarse < 0.0)
  {
    ABORT( status,
         COARSEGRAINFREQUENCYSERIESH_ENEGFMIN,
         COARSEGRAINFREQUENCYSERIESH_MSGENEGFMIN );
  }

  if (fMinFine < 0.0)
  {
    ABORT( status,
         COARSEGRAINFREQUENCYSERIESH_ENEGFMIN,
         COARSEGRAINFREQUENCYSERIESH_MSGENEGFMIN );
  }

  /*    frequency spacing must be positive */

  ASSERT(deltaFCoarse > 0.0, status,
         COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
         COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

  ASSERT(deltaFFine > 0.0, status,
         COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
         COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

  /* check for length mismatch */

  if (output->data->length != lengthCoarse)
  {
    ABORT(status,
         COARSEGRAINFREQUENCYSERIESH_EMMLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEMMLEN);
  }

  /* Calculate coarse-graining parameters */

  offset = ( f0Coarse - f0Fine ) / deltaFFine;

  resRatio = deltaFCoarse / deltaFFine;

  /* Check that coarse-graining makes sense */

  /* make sure coarse-grained series is not finer than fine-grained */
  if ( resRatio < 1.0 )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /* make sure minimum frequency in coarse-grained series is not
     less than minimum frequency in fine-grained series */
  if ( fMinCoarse < fMinFine )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /* make sure maximum frequency in coarse-grained series is not
     more than maximum frequency in fine-grained series */
  if ( offset + resRatio * ( (REAL8) lengthCoarse - 0.5 )
       > lengthFine - 0.5 )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  if (f0Coarse == 0.0)
  {
    /* DC component */

    lamMax = (resRatio / 2.0) - 0.5 ;
    lMax = (UINT4) floor(lamMax);

    if ( lamMax != (REAL8) lMax )
    {
      value = ( lamMax - (REAL8) lMax ) * input->data->data[lMax+1];
    }
    else {
      value = 0.0;
    }

    for ( l = 1 ; l <= lMax ; ++l)
    {
      value += input->data->data[l];
    }

    output->data->data[0] = ( input->data->data[0] + 2.0 * value )
      / resRatio;

    k = 1;
  }
  else
  {
    k = 0;
  }

  for ( ; k < lengthCoarse ; ++k )
  {
    lamMin = offset + ( (REAL8) k - 0.5 ) * resRatio + 0.5 ;
    lMin = (UINT4) ceil(lamMin);

    if ( lamMin != (REAL8) lMin ) {
      value = ( (REAL8) lMin - lamMin ) * input->data->data[lMin-1];
    }
    else
    {
      value = 0.0;
    }

    lamMax = offset + ( (REAL8) k + 0.5 ) * resRatio - 0.5 ;
    lMax = (UINT4) floor(lamMax);

    for ( l = lMin ; l <= lMax ; ++l)
    {
      value += input->data->data[l];
    }

    if ( lamMax != (REAL8) lMax ) {
      value += ( lamMax - (REAL8) lMax ) * input->data->data[lMax+1];
    }

    output->data->data[k] = value / resRatio;

  }

  /* Set output properties */
  output->sampleUnits = input->sampleUnits;
  strncpy( output->name, input->name, LALNameLength );
  output->f0 = f0Coarse;
  output->deltaF = deltaFCoarse;
  output->epoch = input->epoch;

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* LALSCoarseGrainFrequencySeries() */


void
LALDCoarseGrainFrequencySeries(LALStatus                      *status,
                               REAL8FrequencySeries           *output,
                               const REAL8FrequencySeries     *input,
                               const FrequencySamplingParams  *params)

{
  UINT4         lengthCoarse, lengthFine;
  REAL8         f0Coarse, f0Fine;
  REAL8         fMinCoarse, fMinFine;
  REAL8         deltaFCoarse, deltaFFine;
  REAL8         offset, resRatio;

  UINT4         k, l;
  UINT4         lMin, lMax;
  REAL8         lamMin, lamMax;
  REAL8         value;

  /* initialize status structure */
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  /* checks for null pointers: */

  /*    output series */
  ASSERT(output != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of output series */
  ASSERT(output->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of data member of output series */
  ASSERT(output->data->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    input series */
  ASSERT(input != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of input series */
  ASSERT(input->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of data member of input series */
  ASSERT(input->data->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    parameter structure */
  ASSERT(params != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /* checks for duplicate pointers: */
  ASSERT(output != input, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  ASSERT(output->data != input->data, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  ASSERT(output->data->data != input->data->data, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  /* extract coarse-grained parameters  */
  lengthCoarse = params->length;
  f0Coarse     = params->f0;
  deltaFCoarse = params->deltaF;
  fMinCoarse = ( f0Coarse == 0.0 ? 0.0 : f0Coarse - deltaFCoarse / 2.0 );

  /* extract fine-grained parameters */
  lengthFine   = input->data->length;
  f0Fine       = input->f0;
  deltaFFine   = input->deltaF;
  fMinFine = ( f0Fine == 0.0 ? 0.0 : f0Fine - deltaFFine / 2.0 );

  /* check for legality of values */

  /*    length must be positive */

  ASSERT(lengthCoarse != 0, status,
         COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

  ASSERT(lengthFine != 0, status,
         COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

  /*    start frequency must not be negative */

  if (fMinCoarse < 0.0)
  {
    ABORT( status,
         COARSEGRAINFREQUENCYSERIESH_ENEGFMIN,
         COARSEGRAINFREQUENCYSERIESH_MSGENEGFMIN );
  }

  if (fMinFine < 0.0)
  {
    ABORT( status,
         COARSEGRAINFREQUENCYSERIESH_ENEGFMIN,
         COARSEGRAINFREQUENCYSERIESH_MSGENEGFMIN );
  }

  /*    frequency spacing must be positive */

  ASSERT(deltaFCoarse > 0.0, status,
         COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
         COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

  ASSERT(deltaFFine > 0.0, status,
         COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
         COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

  /* check for length mismatch */

  if (output->data->length != lengthCoarse)
  {
    ABORT(status,
         COARSEGRAINFREQUENCYSERIESH_EMMLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEMMLEN);
  }

  /* Calculate coarse-graining parameters */

  offset = ( f0Coarse - f0Fine ) / deltaFFine;

  resRatio = deltaFCoarse / deltaFFine;

  /* Check that coarse-graining makes sense */

  /* make sure coarse-grained series is not finer than fine-grained */
  if ( resRatio < 1.0 )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /* make sure minimum frequency in coarse-grained series is not
     less than minimum frequency in fine-grained series */
  if ( fMinCoarse < fMinFine )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /* make sure maximum frequency in coarse-grained series is not
     more than maximum frequency in fine-grained series */
  if ( offset + resRatio * ( (REAL8) lengthCoarse - 0.5 )
       > lengthFine - 0.5 )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  if (f0Coarse == 0.0)
  {
    /* DC component */

    lamMax = (resRatio / 2.0) - 0.5 ;
    lMax = (UINT4) floor(lamMax);

    if ( lamMax != (REAL8) lMax )
    {
      value = ( lamMax - (REAL8) lMax ) * input->data->data[lMax+1];
    }
    else {
      value = 0.0;
    }

    for ( l = 1 ; l <= lMax ; ++l)
    {
      value += input->data->data[l];
    }

    output->data->data[0] = ( input->data->data[0] + 2.0 * value )
      / resRatio;

    k = 1;
  }
  else
  {
    k = 0;
  }

  for ( ; k < lengthCoarse ; ++k )
  {
    lamMin = offset + ( (REAL8) k - 0.5 ) * resRatio + 0.5 ;
    lMin = (UINT4) ceil(lamMin);

    if ( lamMin != (REAL8) lMin ) {
      value = ( (REAL8) lMin - lamMin ) * input->data->data[lMin-1];
    }
    else
    {
      value = 0.0;
    }

    lamMax = offset + ( (REAL8) k + 0.5 ) * resRatio - 0.5 ;
    lMax = (UINT4) floor(lamMax);

    for ( l = lMin ; l <= lMax ; ++l)
    {
      value += input->data->data[l];
    }

    if ( lamMax != (REAL8) lMax ) {
      value += ( lamMax - (REAL8) lMax ) * input->data->data[lMax+1];
    }

    output->data->data[k] = value / resRatio;

  }

  /* Set output properties */
  output->sampleUnits = input->sampleUnits;
  strncpy( output->name, input->name, LALNameLength );
  output->f0 = f0Coarse;
  output->deltaF = deltaFCoarse;
  output->epoch = input->epoch;

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* LALDCoarseGrainFrequencySeries() */


void
LALCCoarseGrainFrequencySeries(LALStatus                      *status,
                               COMPLEX8FrequencySeries        *output,
                               const COMPLEX8FrequencySeries  *input,
                               const FrequencySamplingParams  *params)

{
  UINT4         lengthCoarse, lengthFine;
  REAL8         f0Coarse, f0Fine;
  REAL8         fMinCoarse, fMinFine;
  REAL8         deltaFCoarse, deltaFFine;
  REAL8         offset, resRatio;

  UINT4         k, l;
  UINT4         lMin, lMax;
  REAL8         lamMin, lamMax;
  COMPLEX8      value;

  /* initialize status structure */
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  /* printf("entering function\n"); */

  /* checks for null pointers: */

  /*    output series */
  ASSERT(output != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of output series */
  ASSERT(output->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of data member of output series */
  ASSERT(output->data->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    input series */
  ASSERT(input != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of input series */
  ASSERT(input->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of data member of input series */
  ASSERT(input->data->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    parameter structure */
  ASSERT(params != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /* checks for duplicate pointers: */

  ASSERT(output != input, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  ASSERT(output->data != input->data, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  ASSERT(output->data->data != input->data->data, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  /* extract coarse-grained parameters  */
  lengthCoarse = params->length;
  f0Coarse     = params->f0;
  deltaFCoarse = params->deltaF;
  fMinCoarse = ( f0Coarse == 0.0 ? 0.0 : f0Coarse - deltaFCoarse / 2.0 );

  /* extract fine-grained parameters */
  lengthFine   = input->data->length;
  f0Fine       = input->f0;
  deltaFFine   = input->deltaF;
  fMinFine = ( f0Fine == 0.0 ? 0.0 : f0Fine - deltaFFine / 2.0 );

  /* check for legality of values */

  /*    length must be positive */

  ASSERT(lengthCoarse != 0, status,
         COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

  ASSERT(lengthFine != 0, status,
         COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

  /*    start frequency must not be negative */

  if (fMinCoarse < 0.0)
  {
    ABORT( status,
         COARSEGRAINFREQUENCYSERIESH_ENEGFMIN,
         COARSEGRAINFREQUENCYSERIESH_MSGENEGFMIN );
  }

  if (fMinFine < 0.0)
  {
    ABORT( status,
         COARSEGRAINFREQUENCYSERIESH_ENEGFMIN,
         COARSEGRAINFREQUENCYSERIESH_MSGENEGFMIN );
  }

  /*    frequency spacing must be positive */

  ASSERT(deltaFCoarse > 0.0, status,
         COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
         COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

  ASSERT(deltaFFine > 0.0, status,
         COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
         COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

  /* check for length mismatch */

  if (output->data->length != lengthCoarse)
  {
    ABORT(status,
         COARSEGRAINFREQUENCYSERIESH_EMMLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEMMLEN);
  }

  /* Calculate coarse-graining parameters */

  offset = ( f0Coarse - f0Fine ) / deltaFFine;

  resRatio = deltaFCoarse /deltaFFine;

  /* Check that coarse-graining makes sense */

  /* make sure coarse-grained series is not finer than fine-grained */
  if ( resRatio < 1.0 )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /* make sure minimum frequency in coarse-grained series is not
     less than minimum frequency in fine-grained series */
  if ( fMinCoarse < fMinFine )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /* make sure maximum frequency in coarse-grained series is not
     more than maximum frequency in fine-grained series */
  if ( offset + resRatio * ( (REAL8) lengthCoarse - 0.5 )
       > lengthFine - 0.5 )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /*
  printf("survived checks\n");

  printf("res ratio %f/%f = %f\n",deltaFCoarse,deltaFFine,resRatio);
  printf("offset (%f-%f)/%f = %f\n",f0Coarse,f0Fine,deltaFFine,offset);
  */

  if (f0Coarse == 0.0)
  {
    /* DC component */

    lamMax = (resRatio / 2.0) - 0.5 ;
    lMax = (UINT4) floor(lamMax);

    /* printf("%f %d %d\n",lamMax,lMax,lengthFine); */

    if ( lamMax != (REAL8) lMax )
    {
      value = ( lamMax - (REAL8) lMax ) * creal(input->data->data[lMax+1]);
    }
    else {
      value = 0.0;
    }

    for ( l = 1 ; l <= lMax ; ++l)
    {
      value += creal(input->data->data[l]);
    }

    output->data->data[0] = creal( input->data->data[0] + 2.0 * value )
      / resRatio;

    /* :TODO: ?  check that imaginary parts of DC vanish? */

    k = 1;
  }
  else
  {
    k = 0;
  }

  for ( ; k < lengthCoarse ; ++k )
  {
    lamMin = offset + ( (REAL8) k - 0.5 ) * resRatio + 0.5 ;
    lMin = (UINT4) ceil(lamMin);

    /* printf("%f %d\n",lamMin,lMin); */

    if ( lamMin != (REAL8) lMin ) {
      value = ( (REAL8) lMin - lamMin ) * input->data->data[lMin-1];
    }
    else
    {
      value = 0.0;
    }

    lamMax = offset + ( (REAL8) k + 0.5 ) * resRatio - 0.5 ;
    lMax = (UINT4) floor(lamMax);

    /*    printf("%f %d %d\n",lamMax,lMax,lengthFine); */

    for ( l = lMin ; l <= lMax ; ++l)
    {
      value += input->data->data[l];
    }

    if ( lamMax != (REAL8) lMax ) {
      value += ( lamMax - (REAL8) lMax ) * input->data->data[lMax+1];
    }

    output->data->data[k] = value / resRatio;

  }

  /* Set output properties */
  output->sampleUnits = input->sampleUnits;
  strncpy( output->name, input->name, LALNameLength );
  output->f0 = f0Coarse;
  output->deltaF = deltaFCoarse;
  output->epoch = input->epoch;

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* LALCCoarseGrainFrequencySeries() */


void
LALZCoarseGrainFrequencySeries(LALStatus                      *status,
                               COMPLEX16FrequencySeries        *output,
                               const COMPLEX16FrequencySeries  *input,
                               const FrequencySamplingParams  *params)

{
  UINT4         lengthCoarse, lengthFine;
  REAL8         f0Coarse, f0Fine;
  REAL8         fMinCoarse, fMinFine;
  REAL8         deltaFCoarse, deltaFFine;
  REAL8         offset, resRatio;

  UINT4         k, l;
  UINT4         lMin, lMax;
  REAL8         lamMin, lamMax;
  COMPLEX16      value;

  /* initialize status structure */
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  /* printf("entering function\n"); */

  /* checks for null pointers: */

  /*    output series */
  ASSERT(output != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of output series */
  ASSERT(output->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of data member of output series */
  ASSERT(output->data->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    input series */
  ASSERT(input != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of input series */
  ASSERT(input->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of data member of input series */
  ASSERT(input->data->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    parameter structure */
  ASSERT(params != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /* checks for duplicate pointers: */

  ASSERT(output != input, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  ASSERT(output->data != input->data, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  ASSERT(output->data->data != input->data->data, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  /* extract coarse-grained parameters  */
  lengthCoarse = params->length;
  f0Coarse     = params->f0;
  deltaFCoarse = params->deltaF;
  fMinCoarse = ( f0Coarse == 0.0 ? 0.0 : f0Coarse - deltaFCoarse / 2.0 );

  /* extract fine-grained parameters */
  lengthFine   = input->data->length;
  f0Fine       = input->f0;
  deltaFFine   = input->deltaF;
  fMinFine = ( f0Fine == 0.0 ? 0.0 : f0Fine - deltaFFine / 2.0 );

  /* check for legality of values */

  /*    length must be positive */

  ASSERT(lengthCoarse != 0, status,
         COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

  ASSERT(lengthFine != 0, status,
         COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

  /*    start frequency must not be negative */

  if (fMinCoarse < 0.0)
  {
    ABORT( status,
         COARSEGRAINFREQUENCYSERIESH_ENEGFMIN,
         COARSEGRAINFREQUENCYSERIESH_MSGENEGFMIN );
  }

  if (fMinFine < 0.0)
  {
    ABORT( status,
         COARSEGRAINFREQUENCYSERIESH_ENEGFMIN,
         COARSEGRAINFREQUENCYSERIESH_MSGENEGFMIN );
  }

  /*    frequency spacing must be positive */

  ASSERT(deltaFCoarse > 0.0, status,
         COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
         COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

  ASSERT(deltaFFine > 0.0, status,
         COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
         COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

  /* check for length mismatch */

  if (output->data->length != lengthCoarse)
  {
    ABORT(status,
         COARSEGRAINFREQUENCYSERIESH_EMMLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEMMLEN);
  }

  /* Calculate coarse-graining parameters */

  offset = ( f0Coarse - f0Fine ) / deltaFFine;

  resRatio = deltaFCoarse /deltaFFine;

  /* Check that coarse-graining makes sense */

  /* make sure coarse-grained series is not finer than fine-grained */
  if ( resRatio < 1.0 )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /* make sure minimum frequency in coarse-grained series is not
     less than minimum frequency in fine-grained series */
  if ( fMinCoarse < fMinFine )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /* make sure maximum frequency in coarse-grained series is not
     more than maximum frequency in fine-grained series */
  if ( offset + resRatio * ( (REAL8) lengthCoarse - 0.5 )
       > lengthFine - 0.5 )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /*
  printf("survived checks\n");

  printf("res ratio %f/%f = %f\n",deltaFCoarse,deltaFFine,resRatio);
  printf("offset (%f-%f)/%f = %f\n",f0Coarse,f0Fine,deltaFFine,offset);
  */

  if (f0Coarse == 0.0)
  {
    /* DC component */

    lamMax = (resRatio / 2.0) - 0.5 ;
    lMax = (UINT4) floor(lamMax);

    /* printf("%f %d %d\n",lamMax,lMax,lengthFine); */

    if ( lamMax != (REAL8) lMax )
    {
      value = creal(( lamMax - (REAL8) lMax ) * input->data->data[lMax+1]);
    }
    else {
      value = 0.0;
    }

    for ( l = 1 ; l <= lMax ; ++l)
    {
      value += creal(input->data->data[l]);
    }

    output->data->data[0] = creal( input->data->data[0] + 2.0 * value )
      / resRatio;

    /* :TODO: ?  check that imaginary parts of DC vanish? */

    k = 1;
  }
  else
  {
    k = 0;
  }

  for ( ; k < lengthCoarse ; ++k )
  {
    lamMin = offset + ( (REAL8) k - 0.5 ) * resRatio + 0.5 ;
    lMin = (UINT4) ceil(lamMin);

    /* printf("%f %d\n",lamMin,lMin); */

    if ( lamMin != (REAL8) lMin ) {
      value = ( (REAL8) lMin - lamMin ) * input->data->data[lMin-1];
    }
    else
    {
      value = 0.0;
    }

    lamMax = offset + ( (REAL8) k + 0.5 ) * resRatio - 0.5 ;
    lMax = (UINT4) floor(lamMax);

    /*    printf("%f %d %d\n",lamMax,lMax,lengthFine); */

    for ( l = lMin ; l <= lMax ; ++l)
    {
      value += input->data->data[l];
    }

    if ( lamMax != (REAL8) lMax ) {
      value += ( lamMax - (REAL8) lMax ) * input->data->data[lMax+1];
    }

    output->data->data[k] = value / resRatio;

  }

  /* Set output properties */
  output->sampleUnits = input->sampleUnits;
  strncpy( output->name, input->name, LALNameLength );
  output->f0 = f0Coarse;
  output->deltaF = deltaFCoarse;
  output->epoch = input->epoch;

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* LALZCoarseGrainFrequencySeries() */
