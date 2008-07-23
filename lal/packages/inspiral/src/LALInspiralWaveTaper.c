/*
*  Copyright (C) 2007 David McKechan, Thomas Cokelaer
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

/*  <lalVerbatim file="LALInspiralWaveTaperCV">
 *  Author:  McKechan D J A
 *  </lalVerbatim>  */

/*  <lalLaTeX>
 *
 *  \subsection{Module \texttt{LALInspiralWaveTaper.c} }
 *
 *  The code \texttt{LALInspiralWaveTaper} imposes a smooth time tapering at
 *  the beginning and/or the end of waves in the time domain. 
 *
 *  It takes in a Real4Vector and searches for the beginning and end points 
 *  of the signal, in case there are null data points at either end. It then
 *  tapers the wave from the ends to the second maxima in the waveform, 
 *  according to formula 3.35 of gr-qc/0001023.
 *
 *  If the waveform does has less than 4 peaks such that it cannot be tapered
 *  from the both ends to the second peak then the central point is tapered to 
 *  insted.
 *
 *  The bookends option allows the user to specify whether just the start, 
 *  just the end or both the start and end of the signal are tapered. 
 *  This corresponds to bookends = 1, 2 or 3 respectively. If bookends does 
 *  not equal 1, 2 or 3 then option 3 is assumed.
 *
 *  \subsubsection*{Prototypes}
 *  \vspace{0.1in}
 *  \input{LALInspiralWaveTaperCP}
 *  \index{\verb&LALInspiralWaveTaper()&}
 *
 *  \subsubsection*{Description}
 *
 *
 *  \subsubsection*{Uses}
 *
 *   
 *
 *   \subsubsection*{Notes}
 *
 *   \vfill{\footnotesize\input{LALInspiralWaveTaperCV}}
 *
 *   </lalLaTeX>  */

#include <lal/LALInspiral.h>
#include <stdio.h>
#include <math.h>

NRCSID(LALINSPIRALWAVETAPERC, 
		  "$Id$"); 

/*  <lalVerbatim file="LALInspiralWaveTaperCP"> */
void LALInspiralWaveTaper(
                   LALStatus   *status,
                   REAL4Vector *signal,
                   UINT4       bookends)
{ /* </lalVerbatim>  */
 
  UINT4 i, start, end, mid, n; /* indices */
  UINT4 flag, safe = 1;
  UINT4 length;   
  REAL4 z, sigma;          
  REAL4 realN, realI;  /* REAL4 values of n & i used for the calculations */

  INITSTATUS(status, "LALInspiralWaveTaper",LALINSPIRALWAVETAPERC);
  ATTATCHSTATUSPTR(status);

  ASSERT(signal, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(signal->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL); 

  length = signal->length;

  if(bookends < 1 || bookends > 3)
    bookends = 3;    

  
  /* Search for start and end of signal */
  flag = 0;
  i = 0;
  while(flag == 0)
  {
    if( signal->data[i] != 0.)
    {
      start = i;
      flag = 1;
    }
    i++;
  }  
  flag = 0;
  i = length - 1;
  while(flag == 0)
  {
    if( signal->data[i] != 0.)
    {
      end = i;
      flag = 1;
    }
    i--;
  }        
  

  /* Check we have more than 2 data points */
  if((end - start) <= 1)
  {
    LALWarning( status, "Data less than 3 points, cannot taper!" );
    safe = 0;
  }

  if( safe == 1)
  {
    /* Calculate middle point in case of short waveform */    
    mid = (start+end)/2;

    /* If requested search for second peak from start and taper */
    if( bookends != 2 )
    {
      flag = 0;
      i = start+1;
      while( flag < 2 && i != mid)
      {
        if( fabs(signal->data[i]) >= fabs(signal->data[i-1]) )
          if( fabs(signal->data[i]) >= fabs(signal->data[i+1]) )
          {
            if( fabs(signal->data[i]) == fabs(signal->data[i+1]) )
              i++;
            flag++;
            n = i - start;
          }
        i++;
      }
      /* Have we reached the middle? */
      if( flag < 2)
      {
        n = mid;
        /* Was it an even length vector? */
        if( mid*2 < end )
          mid++;
      }

      /* Taper to that point */
      realN = (REAL4)(n);
      signal->data[start] = 0.0;
      for(i=start+1; i < n-1; i++)
      {
        realI = (REAL4)(i);
        z = (realN - 1.0)/realI + (realN - 1.0)/(realI - (realN - 1.0));
        sigma = 1.0/(exp(z) + 1.0);
        signal->data[i] = signal->data[i]*sigma;
      }
    }     
  
    /* If requested search for second peak from end */
    if( bookends > 1 )
    {    
      i = end - 1;
      flag = 0;
      while( flag < 2 && i != mid )
      {
        if( fabs(signal->data[i]) >= fabs(signal->data[i+1]) )
          if( fabs(signal->data[i]) >= fabs(signal->data[i-1]) )
          {
            if( fabs(signal->data[i]) == fabs(signal->data[i-1]) )
              i--;
            flag++;
            n = end - i;
          }
        i--;
      }     
      /* Have we reached the middle? */
      if( flag < 2)
      {
        n = mid;
        /* Was it an even length vector? */
        if( mid*2 < end )
          mid++;
      }

      /* Taper to that point */
      realN = (REAL4)(n);
      signal->data[end] = 0.0;   
      for(i=end-1; i > end-n+1; i--)
      {
        realI = (REAL4)(end - i);
        z = (realN - 1.0)/realI + (realN - 1.0)/(realI - (realN-1.0));
        sigma = 1.0/(exp(z) + 1.0);
        signal->data[i] = signal->data[i]*sigma;
      }
    }
  }

  CHECKSTATUSPTR(status);
  DETATCHSTATUSPTR(status);
  RETURN (status);
}


