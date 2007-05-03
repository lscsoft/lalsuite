/*  <lalVerbatim file="LALInspiralWaveTaperCV">
 *  Author:  McKechan D J A
</lalVerbatim>  */

/*  <lalLaTeX>
 *
 *  \subsection{Module \texttt{LALInspiralWaveTaper.c} }
 *
 *  The code \texttt{LALInspiralWaveTaper} imposes a smooth time window on a chirp for 
 *  waves in the time-domain. 
 *
 *  It takes in a Real4Vector and searches for the beginning and end points of the signal.
 *  in case there are null data points at either end.
 *
 *  It tapers the wave over n data points according to formula 3.35 of the Physical Review
 *  D 63 084036
 *
 *  If the user selects a wave to be tapered over n data points then the 1st and nth data points
 *  are 0 and 1. For instance if n = 5, then there are 3 points for which the data is scaled by the
 *  formula.
 *
 *  The bookends option allows the user to specify whether just the start, just the end
 *  or both the start and end of the signal are tapered. This corresponds to bookends = 1,
 *  2 or 3 respectively. If bookends does not equal 1, 2 or 3 then option 3 is assumed.
 *
 *  \subsubsection*{Prototypes}
 *  \vspace{0.1in}
 *  \input{LALInspiralWaveTaperCP}
 *  \index{\verb&LALInspiralWaveTaper()&}
 *  \begin{itemize}
 *  \end{itemize}
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

NRCSID(LALINSPIRALWAVETAPERC, "$Id$"); 

/*  <lalVerbatim file="LALInspiralWaveTaperCP"> */
void LALInspiralWaveTaper(
	 	LALStatus        *status,
	      	REAL4Vector      *signal,
		UINT4		 n,
		UINT4		 bookends
			    )
{ /* </lalVerbatim>  */

  UINT4 i, start, end; /* indices */
  UINT4 flag, abort;      /* flags */
  UINT4 length;        
  REAL4 z, sigma;          
  REAL4 rn, ri;           /* REAL4 values of n & i used for the calculations */

  INITSTATUS(status, "LALInspiralWaveTaper",LALINSPIRALWAVETAPERC);
  ATTATCHSTATUSPTR(status);

  ASSERT(signal, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(signal->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL); 

  length = signal->length;
 
  rn = (REAL4)(n);

  /* Check options */
  /* minimum value is 2 - recommended is (?)*/
  if(rn < 2)
    rn = 2;
  
  if(bookends < 1 || bookends > 3)
    bookends = 3;	  

  /* Check to see if we can taper the signal. 
   * There maybe null data points at the start or the end of the vector 
   * so we find the start and end points of the index and check that they are n points
   * from the start and end of the vector */ 
  abort = 0;
  
  /* If required search index for start of signal */
  if(bookends != 2)
    flag = 0;
    i = 0;
    while(flag == 0){
      if( signal->data[i] != 0){
        start = i;
        flag = 1;
      }
      i++;
      /* safety */
      if( i > (length - n) ){
        flag = 1;	     
        abort = 1;
      }       
    }     

  /* If required search index for end of signal */
  if( bookends > 1){
    flag = 0;
    i = length - 1;
    while(flag == 0){
      if( signal->data[i] != 0){
        end = i;
        flag = 1;
      }
      i--;
      /* safety */
      if( i < n ){
        flag = 1;	     
        abort = 1;
      }        
    }
  }    

  /* If we can taper the signal then we do it here */ 	   
  /* Otherwise the signal is returned with no changes */
  if( abort == 0){
    /* end points */
    if(bookends != 2)
      signal->data[start] = 0.0; 
    if(bookends > 1)
      signal->data[end] = 0.0;

    for( i = 1; i < (n-1); i++){
      ri = (REAL4)(i);
      
      z = (rn - 1.0)/(ri - (REAL4)(start)) + (rn - 1.0)/(ri - (rn-1.0));  
      sigma = 1.0/(exp(z) + 1.0);
      
      if(bookends != 2)
        signal->data[start+i] *= sigma;
      if(bookends > 1)
        signal->data[end-i] *= sigma;
    }     
    
  }     


  CHECKSTATUSPTR(status);
	    
  DETATCHSTATUSPTR(status);
  RETURN (status);
}


