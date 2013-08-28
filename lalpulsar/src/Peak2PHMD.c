/*
 *  Copyright (C) 2005 Badri Krishnan, Alicia Sintes
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

/*
 *
 * History:   Created by Sintes June 20, 2001
 *            Modified...
 *
 *
*/

#include <lal/PHMD.h>


/**
 * \author Sintes, A. M.
 * \brief Construction of Partial-Hough-Map-Derivatives (\c phmd) given a peak-gram and the look-up-table.
 * \ingroup PHMD_h
 *
 * ### Description ###
 *
 * This routine produces a \c phmd at a certain frequency for a given  peak-gram and
 * look-up-table.
 *
 * The inputs are:
 *
 * <tt>phmd->fBin</tt>: The frequency bin of this \c phmd.
 *
 * <tt>*lut</tt>: The look-up-table (of type  \c HOUGHptfLUT)
 *
 * <tt>*pg</tt>: The peak-gram  (of type  \c HOUGHPeakGram)
 *
 * The function LALHOUGHPeak2PHMD() makes sure that the  \c lut, the
 * peak-gram and also the frequency of the \c phmd
 * are compatible.
 *
 * The output <tt>HOUGHphmd  *phmd</tt> is  a structure
 * containing the frequency bin of this \c phmd,
 * the total number of borders of each type (<em>Left and Right</em>) to be
 * marked, the pointers to the borders in the corresponding
 * look-up-table, plus  \e border effects of clipping  on a finite
 * patch.
 *
 */
void LALHOUGHPeak2PHMD (LALStatus    *status,
			HOUGHphmd    *phmd, /* partial Hough map derivative */
			HOUGHptfLUT  *lut, /* Look up table */
			HOUGHPeakGram *pg)  /* peakgram */
{

  INT2    i,j;
  UCHAR   *column1P;
  INT4    fBinDif,shiftPeak,minPeakBin,maxPeakBin,thisPeak;
  INT2    relatIndex;
  UINT4   lengthLeft,lengthRight,n, searchIndex;
  UINT8   firstBin,lastBin,pgI,pgF;
  /* --------------------------------------------- */

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

/**
 * lots of error checking of arguments -- asserts have been
 * replaced by aborts
 */

  /*   Make sure the arguments are not NULL: */
  if (phmd == NULL) {
    /* fprintf(stderr,"null pointer found [Peak2PHMD.c %d]\n", __LINE__); */
    ABORT( status, PHMDH_ENULL, PHMDH_MSGENULL);
  }
  /* ASSERT (phmd, status, PHMDH_ENULL, PHMDH_MSGENULL); */

  if (lut == NULL) {
    /* fprintf(stderr,"null pointer found [Peak2PHMD.c %d]\n", __LINE__); */
    ABORT( status, PHMDH_ENULL, PHMDH_MSGENULL);
  }
  /* ASSERT (lut,  status, PHMDH_ENULL, PHMDH_MSGENULL); */

  if (pg == NULL) {
    /* fprintf(stderr,"null pointer found [Peak2PHMD.c %d]\n", __LINE__); */
    ABORT( status, PHMDH_ENULL, PHMDH_MSGENULL);
  }
  /* ASSERT (pg,   status, PHMDH_ENULL, PHMDH_MSGENULL); */

  if (phmd->firstColumn == NULL) {
    /* fprintf(stderr,"null pointer found [Peak2PHMD.c %d]\n", __LINE__); */
    ABORT( status, PHMDH_ENULL, PHMDH_MSGENULL);
  }
  /* ASSERT (phmd->firstColumn, status, PHMDH_ENULL, PHMDH_MSGENULL); */

  /*  Make sure there are elements in firstColumn */
  if (phmd->ySide == 0) {
    /* fprintf(stderr,"phmd->ySide should be non-zero [Peak2PHMD.c %d]\n", __LINE__); */
    ABORT( status, PHMDH_ESIZE, PHMDH_MSGESIZE);
  }
  /* ASSERT (phmd->ySide, status, PHMDH_ESIZE, PHMDH_MSGESIZE); */

  /* Make sure peakgram and lut have same frequency discretization */
  if ( fabs((REAL4)lut->deltaF - (REAL4)pg->deltaF) > 1.0e-6) {
    ABORT( status, PHMDH_EVAL, PHMDH_MSGEVAL);
  }
  /* ASSERT ( fabs((REAL4)lut->deltaF - (REAL4)pg->deltaF) < 1.0e-6,  status, PHMDH_EVAL, PHMDH_MSGEVAL); */

  /* Make sure phmd.fBin and lut are compatible */
  fBinDif = abs( (phmd->fBin) - (lut->f0Bin) );
  if ( fBinDif > lut->nFreqValid ) {
    /* fprintf(stderr,"fBinDif > nFreqValid [Peak2PHMD.c %d]\n", __LINE__); */
    ABORT( status, PHMDH_EFREQ, PHMDH_MSGEFREQ);
  }
  /* ASSERT (fBinDif < lut->nFreqValid, status, PHMDH_EFREQ, PHMDH_MSGEFREQ); */


  pgI = pg->fBinIni;
  pgF = pg->fBinFin;
  /* bounds of interval to look at in the peakgram */
  firstBin = (phmd->fBin) + (lut->iniBin) + (lut->offset);
  lastBin  = firstBin + (lut->nBin)-1;

  /* Make sure peakgram f-interval and phmd.fBin+lut are compatible */
  /*   ASSERT ( pgI <= firstBin, status, PHMDH_EINT, PHMDH_MSGEINT); */
  /*   ASSERT ( pgF >= lastBin,  status, PHMDH_EINT, PHMDH_MSGEINT); */
  if (  pgI > firstBin ) {
    ABORT( status, PHMDH_EINT, PHMDH_MSGEINT);
  }
  if (  pgF < lastBin ) {
    ABORT( status, PHMDH_EINT, PHMDH_MSGEINT);
  }


  /* -------------------------------------------------------------------   */
  /* initialization */
  /* -------------------------------------------------------------------   */
  n= pg->length;

  lengthLeft = 0;
  lengthRight= 0;

  column1P = &(phmd->firstColumn[0]);

  for(i=0; i< phmd->ySide; ++i ){
    *column1P = 0;
    ++column1P;
  }

  minPeakBin = firstBin - pgI;
  maxPeakBin =  lastBin - pgI;

 /* -------------------------------------------------------------------   */
         /* -------------------------------------------   */
  if(n){  /* only if there are peaks present */
         /* -------------------------------------------   */
    INT2   lb1,rb1,lb2,rb2;
    INT2   max1,min1,max2,min2;
    INT2   nBinPos;
    UCHAR  test1;

    nBinPos = (lut->iniBin) + (lut->nBin) -1;
    shiftPeak = pgI - (phmd->fBin) - (lut->offset);

   /* -------------------------------------------------------------------   */
   /* searching for the initial peak to look at */
   /* -------------------------------------------------------------------   */

    searchIndex = (n * minPeakBin ) /(pgF-pgI+1);
    if (searchIndex >= n) { searchIndex = n-1;}

    /* moving backwards */
    /* -----------------*/

    test1=1;
    while (searchIndex >0 && test1) {
      if ( pg->peak[searchIndex -1] >=  minPeakBin ) {
        --searchIndex;
      }
      else {
        test1=0;
      }
    }

    /* moving forwards */
    /* -----------------*/

    test1=1;
    while (searchIndex < n-1  && test1) {
      if ( pg->peak[searchIndex] <  minPeakBin ) {
        ++searchIndex;
      }
      else {
        test1=0;
      }
    }

    /* -------------------------------------------------------------------   */
    /* for all the interesting peaks (or none) */
    /* -------------------------------------------------------------------   */

    test1=1;

    while (searchIndex < n  && test1) {
      thisPeak = pg->peak[searchIndex];

      if (thisPeak > maxPeakBin) {
        test1=0;
      } else {
        if ( thisPeak >= minPeakBin) { /* we got a peak */
	  /* relative Index */
	  relatIndex = thisPeak + shiftPeak;

	  i = relatIndex;
	  if( relatIndex < 0 )  i =  nBinPos - relatIndex;


	  if ( i >= lut->nBin ) {
	    fprintf(stderr,"current index i=%d not lesser than nBin=%d\n [Peak2PHMD.c %d]\n", i,lut->nBin, __LINE__);
	  }


	  /* Reading the bin information */
	  lb1 = lut->bin[i].leftB1;
	  rb1 = lut->bin[i].rightB1;
	  lb2 = lut->bin[i].leftB2;
	  rb2 = lut->bin[i].rightB2;

	  max1 = lut->bin[i].piece1max;
	  min1 = lut->bin[i].piece1min;
	  max2 = lut->bin[i].piece2max;
	  min2 = lut->bin[i].piece2min;

	  /* border selection from lut */

	  if(lb1){
	    if (  lut->border + lb1 == NULL ) {
	      /* fprintf(stderr,"null pointer found [Peak2PHMD.c %d]\n", __LINE__); */
	      ABORT( status, PHMDH_ENULL, PHMDH_MSGENULL);
	    }
	    phmd->leftBorderP[lengthLeft] = &( lut->border[lb1] );
	    ++lengthLeft;
	  }


	  if(lb2){
	    if (  lut->border + lb2 == NULL ) {
	      /* fprintf(stderr,"null pointer found [Peak2PHMD.c %d]\n", __LINE__); */
	      ABORT( status, PHMDH_ENULL, PHMDH_MSGENULL);
	    }
	    phmd->leftBorderP[lengthLeft] = &( lut->border[lb2] );
	    ++lengthLeft;
	  }

	  if(rb1){
	    if (  lut->border + rb1 == NULL ) {
	      /* fprintf(stderr,"null pointer found [Peak2PHMD.c %d]\n", __LINE__); */
	      ABORT( status, PHMDH_ENULL, PHMDH_MSGENULL);
	    }
	    phmd->rightBorderP[lengthRight] = &( lut->border[rb1] );
	    ++lengthRight;
	  }

	  if(rb2){
	    if (  lut->border + rb2 == NULL ) {
	      /* fprintf(stderr,"null pointer found [Peak2PHMD.c %d]\n", __LINE__); */
	      ABORT( status, PHMDH_ENULL, PHMDH_MSGENULL);
	    }
	    phmd->rightBorderP[lengthRight] = &( lut->border[rb2] );
	    ++lengthRight;
	  }

	  /* correcting 1st column */
	  for(j=min1; j<=max1; ++j) {
	    if (phmd->firstColumn + j == NULL) {
	      /* fprintf(stderr,"null pointer found [Peak2PHMD.c %d]\n", __LINE__); */
	      ABORT( status, PHMDH_ENULL, PHMDH_MSGENULL);
	    }
	    phmd->firstColumn[j] = 1;
	  }

	  for(j=min2; j<=max2; ++j) {
	    if (phmd->firstColumn + j == NULL) {
	      /* fprintf(stderr,"null pointer found [Peak2PHMD.c %d]\n", __LINE__); */
	      ABORT( status, PHMDH_ENULL, PHMDH_MSGENULL);
	    }
	    phmd->firstColumn[j] = 1;
	  }

	}
        ++searchIndex;
      }
      /* -------------------------------------------------------------------   */

    }
  }

  /* -------------------------------------------------------------------   */

  phmd->lengthLeft = lengthLeft;
  phmd->lengthRight= lengthRight;
  /* -------------------------------------------   */

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}
