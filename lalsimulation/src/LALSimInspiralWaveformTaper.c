/*
*  Copyright (C) 2007 David McKechan, Thomas Cokelaer, Drew Keppel
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
 * \author  McKechan D J A
 * \file
 * \ingroup LALSimInspiral_h
 *
 * \brief The code \c XLALInspiralREAL4WaveformTaper() and
 * \c XLALInspiralREAL8WaveTaper() impose a smooth time tapering at the
 * beginning and/or the end of REAL4 or REAL8 waveforms in the time domain.
 *
 * They take either a ::REAL4Vector or a ::REAL8Vector and search for the
 * beginning and end points of the signal, in case there are null data points
 * at either end. They then taper the waveform from the ends to the second
 * maxima from each end in the waveform, according to formula 3.35 of
 * <tt>gr-qc/0001023</tt>.
 *
 * If the waveform does has less than 4 maxima, such that it cannot be tapered
 * from the each end to the second peak then the waveform is tapered from the
 * ends to the centre of the instead.
 *
 * The bookends option is an ::LALSimInspiralApplyTaper enumerator and allows the
 * user to specify whether just the start, just the end or both the start and
 * end of the signal are tapered. These options are #LAL_SIM_INSPIRAL_TAPER_START,
 * #LAL_SIM_INSPIRAL_TAPER_END and #LAL_SIM_INSPIRAL_TAPER_STARTEND.
 *
 * \heading{Prototypes}
 * <tt>XLALInspiralREAL4WaveformTaper()</tt>
 * <tt>XLALInspiralREAL8WaveformTaper()</tt>
 *
 * \heading{Description}
 *
 * \heading{Uses}
 *
 * \heading{Notes}
 *
 */

#include <LALSimInspiral.h>
#include <lal/LALError.h>
#include <stdio.h>
#include <math.h>

int XLALSimInspiralREAL4WaveTaper(
		REAL4Vector              *signalvec,	/**< pointer to waveform vector */
		LALSimInspiralApplyTaper  bookends	/**< taper type enumerator */
		)
{
  UINT4 i, start=0, end=0, mid, n=0; /* indices */
  UINT4 flag, safe = 1;
  UINT4 length;
  REAL4 z, sigma;
  REAL4 realN, realI;  /* REAL4 values of n & i used for the calculations */

#ifndef LAL_NDEBUG
  if ( !signalvec )
    XLAL_ERROR( XLAL_EFAULT );

  if ( !signalvec->data )
    XLAL_ERROR( XLAL_EFAULT );
#endif

  /* Check we have chosen a valid tapering method */
  if ( (UINT4) bookends >= (UINT4) LAL_SIM_INSPIRAL_TAPER_NUM_OPTS )
    XLAL_ERROR( XLAL_EINVAL );

  length = signalvec->length;

  if( bookends == LAL_SIM_INSPIRAL_TAPER_NONE )
  {
    XLALPrintWarning( "No taper specified; not tapering.\n" );
    return XLAL_SUCCESS;
  }

  /* Search for start and end of signal */
  flag = 0;
  i = 0;
  while(flag == 0 && i < length )
  {
    if( signalvec->data[i] != 0.)
    {
      start = i;
      flag = 1;
    }
    i++;
  }
  if ( flag == 0 )
  {
    XLALPrintWarning( "No signal found in the vector. Cannot taper.\n" );
    return XLAL_SUCCESS;
  }

  flag = 0;
  i = length - 1;
  while(flag == 0)
  {
    if( signalvec->data[i] != 0.)
    {
      end = i;
      flag = 1;
    }
    i--;
  }


  /* Check we have more than 2 data points */
  if((end - start) <= 1)
  {
    XLALPrintWarning( "Data less than 3 points, cannot taper!\n" );
    safe = 0;
  }

  if( safe == 1)
  {
    /* Calculate middle point in case of short waveform */
    mid = (start+end)/2;

    /* If requested search for second peak from start and taper */
    if( bookends != LAL_SIM_INSPIRAL_TAPER_END )
    {
      flag = 0;
      i = start+1;
      while( flag < 2 && i != mid)
      {
        if( fabs(signalvec->data[i]) >= fabs(signalvec->data[i-1]) )
          if( fabs(signalvec->data[i]) >= fabs(signalvec->data[i+1]) )
          {
            if( fabs(signalvec->data[i]) == fabs(signalvec->data[i+1]) )
              i++;
            flag++;
            n = i - start;
          }
        i++;
      }
      /* Have we reached the middle? */
      if( flag < 2)
        n = mid - start;

      /* Taper to that point */
      realN = (REAL4)(n);
      signalvec->data[start] = 0.0;
      for(i=start+1; i < start + n - 1; i++)
      {
        realI = (REAL4)(i - start);
        z = (realN - 1.0)/realI + (realN - 1.0)/(realI - (realN - 1.0));
        sigma = 1.0/(exp(z) + 1.0);
        signalvec->data[i] = signalvec->data[i]*sigma;
      }
    }

    /* If requested search for second peak from end */
    if( bookends == LAL_SIM_INSPIRAL_TAPER_END || bookends == LAL_SIM_INSPIRAL_TAPER_STARTEND )
    {
      i = end - 1;
      flag = 0;
      while( flag < 2 && i != mid )
      {
        if( fabs(signalvec->data[i]) >= fabs(signalvec->data[i+1]) )
          if( fabs(signalvec->data[i]) >= fabs(signalvec->data[i-1]) )
          {
            if( fabs(signalvec->data[i]) == fabs(signalvec->data[i-1]) )
              i--;
            flag++;
            n = end - i;
          }
        i--;
      }
      /* Have we reached the middle? */
      if( flag < 2)
      {
        n = end - mid;
      }

      /* Taper to that point */
      realN = (REAL4)(n);
      signalvec->data[end] = 0.0;
      for(i=end-1; i > end-n+1; i--)
      {
        realI = (REAL4)(end - i);
        z = (realN - 1.0)/realI + (realN - 1.0)/(realI - (realN-1.0));
        sigma = 1.0/(exp(z) + 1.0);
        signalvec->data[i] = signalvec->data[i]*sigma;
      }
    }
  }

  return XLAL_SUCCESS;
}

int XLALSimInspiralREAL8WaveTaper(
		REAL8Vector              *signalvec,	/**< pointer to waveform vector */
		LALSimInspiralApplyTaper  bookends	/**< taper type enumerator */
		)
{
  UINT4 i, start=0, end=0, mid, n=0; /* indices */
  UINT4 flag, safe = 1;
  UINT4 length;
  REAL8 z, sigma;
  REAL8 realN, realI;  /* REAL4 values of n & i used for the calculations */

#ifndef LAL_NDEBUG
  if ( !signalvec )
    XLAL_ERROR( XLAL_EFAULT );

  if ( !signalvec->data )
    XLAL_ERROR( XLAL_EFAULT );
#endif

  /* Check we have chosen a valid tapering method */
  if ( (UINT4) bookends >= (UINT4) LAL_SIM_INSPIRAL_TAPER_NUM_OPTS )
    XLAL_ERROR( XLAL_EINVAL );

  length = signalvec->length;

  if( bookends == LAL_SIM_INSPIRAL_TAPER_NONE )
  {
    XLALPrintWarning( "No taper specified; not tapering.\n" );
    return XLAL_SUCCESS;
  }

  /* Search for start and end of signal */
  flag = 0;
  i = 0;
  while(flag == 0 && i < length )
  {
    if( signalvec->data[i] != 0.)
    {
      start = i;
      flag = 1;
    }
    i++;
  }
  if ( flag == 0 )
  {
    XLALPrintWarning( "No signal found in the vector. Cannot taper.\n" );
    return XLAL_SUCCESS;
  }

  flag = 0;
  i = length - 1;
  while(flag == 0)
  {
    if( signalvec->data[i] != 0.)
    {
      end = i;
      flag = 1;
    }
    i--;
  }


  /* Check we have more than 2 data points */
  if((end - start) <= 1)
  {
    XLALPrintWarning( "Data less than 3 points, cannot taper!\n" );
    safe = 0;
  }

  if( safe == 1)
  {
    /* Calculate middle point in case of short waveform */
    mid = (start+end)/2;

    /* If requested search for second peak from start and taper */
    if( bookends != LAL_SIM_INSPIRAL_TAPER_END )
    {
      flag = 0;
      i = start+1;
      while( flag < 2 && i != mid)
      {
        if( fabs(signalvec->data[i]) >= fabs(signalvec->data[i-1]) )
          if( fabs(signalvec->data[i]) >= fabs(signalvec->data[i+1]) )
          {
            if( fabs(signalvec->data[i]) == fabs(signalvec->data[i+1]) )
              i++;
            flag++;
            n = i - start;
          }
        i++;
      }
      /* Have we reached the middle? */
      if( flag < 2)
        n = mid - start;

      /* Taper to that point */
      realN = (REAL8)(n);
      signalvec->data[start] = 0.0;
      for(i=start+1; i < start + n - 1; i++)
      {
        realI = (REAL8)(i - start);
        z = (realN - 1.0)/realI + (realN - 1.0)/(realI - (realN - 1.0));
        sigma = 1.0/(exp(z) + 1.0);
        signalvec->data[i] = signalvec->data[i]*sigma;
      }
    }

    /* If requested search for second peak from end */
    if( bookends == LAL_SIM_INSPIRAL_TAPER_END || bookends == LAL_SIM_INSPIRAL_TAPER_STARTEND )
    {
      i = end - 1;
      flag = 0;
      while( flag < 2 && i != mid )
      {
        if( fabs(signalvec->data[i]) >= fabs(signalvec->data[i+1]) )
          if( fabs(signalvec->data[i]) >= fabs(signalvec->data[i-1]) )
          {
            if( fabs(signalvec->data[i]) == fabs(signalvec->data[i-1]) )
              i--;
            flag++;
            n = end - i;
          }
        i--;
      }
      /* Have we reached the middle? */
      if( flag < 2)
      {
        n = end - mid;
      }

      /* Taper to that point */
      realN = (REAL8)(n);
      signalvec->data[end] = 0.0;
      for(i=end-1; i > end-n+1; i--)
      {
        realI = (REAL8)(end - i);
        z = (realN - 1.0)/realI + (realN - 1.0)/(realI - (realN-1.0));
        sigma = 1.0/(exp(z) + 1.0);
        signalvec->data[i] = signalvec->data[i]*sigma;
      }
    }
  }

  return XLAL_SUCCESS;
}
