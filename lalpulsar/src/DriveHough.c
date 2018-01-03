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

#include <lal/LALHough.h>

/** \cond DONT_DOXYGEN */

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/** \endcond */

/*
 * The functions that make up the guts of this module
 */

/** \addtogroup LALHough_h */
/*@{*/

/**
 * constructs the space of <tt>phmd</tt> <tt>PHMDVectorSequence *phmdVS</tt>, given a
 * HOUGHPeakGramVector *pgV and HOUGHptfLUTVector *lutV.
 * The minimum frequency bin present corresponds to <tt>phmdVS->fBinMin</tt> and the
 * total number of different frequencies is
 * <tt>phmdVS->nfSize</tt>. At this moment the \c fBinMin line corresponds to the
 * first row of the cylinder and <tt>phmdVS->breakLine</tt> is set to zero.
 * <tt>phmdVS->breakLine</tt>
 * \f$\in\, [0,\f$\c nfSize) is <em>the pointer</em> which
 * identifies the position of the  \c fBinMin row in the circular-cylinder buffer.
 */
void LALHOUGHConstructSpacePHMD  (LALStatus            *status,	/**< pointer to LALStatus structure */
				  PHMDVectorSequence   *phmdVS, /**< Cylindrical buffer of PHMDs */
				  HOUGHPeakGramVector  *pgV, 	/**< Vetor of peakgrams */
				  HOUGHptfLUTVector    *lutV 	/**< vector of look up tables */)
{

  UINT4    k,j;
  UINT4    nfSize;    /* number of different frequencies */
  UINT4    length;    /* number of elements for each frequency */
  UINT8    fBinMin;   /* present minimum frequency bin */
  UINT8    fBin;      /* present frequency bin */


  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);


  /*   Make sure the arguments are not NULL: */
  ASSERT (phmdVS, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (pgV,    status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (lutV,   status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  ASSERT (phmdVS->phmd, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (pgV->pg,      status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (lutV->lut,    status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  /* Make sure there is no size mismatch */
  ASSERT (pgV->length == lutV->length, status,
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  ASSERT (pgV->length == phmdVS->length, status,
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  /* -------------------------------------------   */

  /* Make sure there are elements to be computed*/
  ASSERT (phmdVS->length, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  ASSERT (phmdVS->nfSize, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);

  /* at the  beggining, the  fBinMin line corresponds to the first row */
  phmdVS->breakLine = 0; /* mark [0,nfSize) (of the circular buffer)
			    pointing to the starting of the fBinMin line */

  length = phmdVS->length;
  nfSize = phmdVS->nfSize;
  fBinMin = phmdVS->fBinMin;
#ifndef LAL_NDEBUG
  REAL8 deltaF =  phmdVS->deltaF = lutV->lut[0].deltaF;   /* frequency resolution */
#endif

  for ( k=0; k<length; ++k ){

    /* make sure all deltaF are consistent */
    ASSERT (deltaF == lutV->lut[k].deltaF,
	    status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);

    fBin = fBinMin;

    for ( j=0; j<  nfSize; ++j ){
      phmdVS->phmd[ j*length+k ].fBin = fBin;

      TRY( LALHOUGHPeak2PHMD(status->statusPtr,
			     &(phmdVS->phmd[ j*length+k ]),
			     &(lutV->lut[k]), &(pgV->pg[k]) ), status);
      ++fBin;
    }
  }


  DETATCHSTATUSPTR (status);
   /* normal exit */
  RETURN (status);
}

/**
 * This function updates the space of <tt>phmd</tt> increasing the frequency <tt>phmdVS->fBinMin</tt> by one.
 */
void LALHOUGHupdateSpacePHMDup  (LALStatus            *status,
				  PHMDVectorSequence   *phmdVS,
				  HOUGHPeakGramVector  *pgV,
				  HOUGHptfLUTVector    *lutV)
{
  UINT4    k,breakLine;
  UINT4    nfSize;    /* number of different frequencies */
  UINT4    length;    /* number of elements for each frequency */
  UINT8    fBinMin;   /* minimum frequency bin */
  UINT8    fBin;      /* present frequency bin */


  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);


  /*   Make sure the arguments are not NULL: */
  ASSERT (phmdVS, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (pgV,    status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (lutV,   status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  ASSERT (phmdVS->phmd, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (pgV->pg,      status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (lutV->lut,    status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  /* Make sure there is no size mismatch */
  ASSERT (pgV->length == lutV->length, status,
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  ASSERT (pgV->length == phmdVS->length, status,
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  /* -------------------------------------------   */

  /* Make sure there are elements to be computed*/
  ASSERT (phmdVS->length, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  ASSERT (phmdVS->nfSize, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
 /* -------------------------------------------   */


  length = phmdVS->length;
  nfSize = phmdVS->nfSize;
#ifndef LAL_NDEBUG
  REAL8 deltaF =  phmdVS->deltaF;
#endif

  breakLine = phmdVS->breakLine; /* old Break Line */
  fBinMin = phmdVS->fBinMin; /* initial frequency value  */

  /* Make sure initial breakLine is in [0,nfSize)  */
  ASSERT ( breakLine < nfSize, status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);

  /* -------------------------------------------   */

  /* Updating the space of PHMD increasing frequency */

  fBin = fBinMin + nfSize;

  for ( k=0; k<length; ++k ){
    /* make sure all deltaF are consistent */
    ASSERT (deltaF == lutV->lut[k].deltaF,
	    status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);

    phmdVS->phmd[ breakLine*length+k ].fBin = fBin;
    TRY( LALHOUGHPeak2PHMD(status->statusPtr,
			   &(phmdVS->phmd[ breakLine*length+k ]),
			   &(lutV->lut[k]), &(pgV->pg[k]) ), status);
  }

  /* Shift fBinMin and its mark */
  ++phmdVS->fBinMin;

  phmdVS->breakLine = (breakLine +1) % nfSize;
  /* mark [0,nfSize) (of the circular buffer, modulus nfSize)
     pointing to the starting of the new fBinMin line */

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}


/**
 * Function for shifting the cylindrical buffer of PHMDs down by one
 * frequency bin -- the highest frequency bin is dropped and an
 * extra frequency bin is added at the lowest frequency
 */
void LALHOUGHupdateSpacePHMDdn  (LALStatus            *status,
				 PHMDVectorSequence   *phmdVS,
				 HOUGHPeakGramVector  *pgV,
				 HOUGHptfLUTVector    *lutV)
{
  UINT4    k,breakLine;
  UINT4    nfSize;    /* number of different frequencies */
  UINT4    length;    /* number of elements for each frequency */

  UINT8    fBin;      /* present frequency bin */

  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);


  /*   Make sure the arguments are not NULL: */
  ASSERT (phmdVS, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (pgV,    status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (lutV,   status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  ASSERT (phmdVS->phmd, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (pgV->pg,      status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (lutV->lut,    status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  /* Make sure there is no size mismatch */
  ASSERT (pgV->length == lutV->length, status,
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  ASSERT (pgV->length == phmdVS->length, status,
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  /* -------------------------------------------   */

  /* Make sure there are elements to be computed*/
  ASSERT (phmdVS->length, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  ASSERT (phmdVS->nfSize, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  /* -------------------------------------------   */

  length = phmdVS->length;
  nfSize = phmdVS->nfSize;
#ifndef LAL_NDEBUG
  REAL8 deltaF =  phmdVS->deltaF;
#endif

  breakLine = phmdVS->breakLine; /* old Break Line */

  /* Make sure initial breakLine is in [0,nfSize)  */
  ASSERT ( breakLine < nfSize, status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);

  /* -------------------------------------------   */

  /* Updating the space of PHMD decreasing frequency */

  /* Shift fBinMin and its mark */
  fBin =  --phmdVS->fBinMin; /* initial frequency value  */

  phmdVS->breakLine = (breakLine + nfSize- 1) % nfSize;
  /* mark [0,nfSize) (of the circular buffer, modulus nfSize)
     pointing to the starting of the new fBinMin line */

  breakLine = phmdVS->breakLine; /* the new Break Line */

  for ( k=0; k<length; ++k ){
    /* make sure all deltaF are consistent */
    ASSERT (deltaF == lutV->lut[k].deltaF,
	    status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);

    phmdVS->phmd[ breakLine*length+k ].fBin = fBin;
    TRY( LALHOUGHPeak2PHMD(status->statusPtr,
			   &(phmdVS->phmd[ breakLine*length+k ]),
			   &(lutV->lut[k]), &(pgV->pg[k]) ), status);
  }



  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}


/**
 * Given PHMDVectorSequence *phmdVS, the space of \c phmd, and
 * UINT8FrequencyIndexVector *freqInd, a structure containing the frequency
 * indices  of the   \c phmd at different time stamps that have to be combined
 * to form a Hough map, the function LALHOUGHConstructHMT() produces the
 * total Hough map.
 */
void LALHOUGHConstructHMT  (LALStatus                  *status,	/**< pointer to LALStatus structure */
			    HOUGHMapTotal              *ht, 	/**< The output hough map */
			    UINT8FrequencyIndexVector  *freqInd,/**< time-frequency trajectory */
			    PHMDVectorSequence         *phmdVS 	/**< set of partial hough map derivatives */)
{


  UINT4    k,j;
  UINT4    breakLine;
  UINT4    nfSize;    /* number of different frequencies */
  UINT4    length;    /* number of elements for each frequency */
  UINT8    fBinMin;   /* present minimum frequency bin */
  INT8     fBin;      /* present frequency bin */
  UINT2    xSide,ySide;

  HOUGHMapDeriv hd; /* the Hough map derivative */

  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (phmdVS,  status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (ht,      status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (freqInd, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  ASSERT (phmdVS->phmd,  status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (freqInd->data, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  /* Make sure there is no size mismatch */
  ASSERT (freqInd->length == phmdVS->length, status,
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  ASSERT (freqInd->deltaF == phmdVS->deltaF, status,
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  /* -------------------------------------------   */

  /* Make sure there are elements  */
  ASSERT (phmdVS->length, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  ASSERT (phmdVS->nfSize, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  /* -------------------------------------------   */

   /* Make sure the ht map contains some pixels */
  ASSERT (ht->xSide, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  ASSERT (ht->ySide, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);

  length = phmdVS->length;
  nfSize = phmdVS->nfSize;

  fBinMin = phmdVS->fBinMin; /* initial frequency value  od the cilinder*/

  breakLine = phmdVS->breakLine;

  /* number of physical pixels */
  xSide = ht->xSide;
  ySide = ht->ySide;

  /* Make sure initial breakLine is in [0,nfSize)  */
  ASSERT ( breakLine < nfSize, status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);

  /* -------------------------------------------   */

  /* Initializing  hd map and memory allocation */
  hd.xSide = xSide;
  hd.ySide = ySide;
  hd.map = (HoughDT *)LALMalloc(ySide*(xSide+1)*sizeof(HoughDT));
  if (hd. map == NULL) {
    ABORT( status, LALHOUGHH_EMEM, LALHOUGHH_MSGEMEM);
  }
  /* -------------------------------------------   */


  TRY( LALHOUGHInitializeHD(status->statusPtr, &hd), status);
  for ( k=0; k<length; ++k ){
    /* read the frequency index and make sure is in the proper interval*/
    fBin =freqInd->data[k] -fBinMin;

    ASSERT ( fBin < nfSize, status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);
    ASSERT ( fBin >= 0,     status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);

    /* find index */
    j = (fBin + breakLine) % nfSize;

    /* Add the corresponding PHMD to HD */
    TRY( LALHOUGHAddPHMD2HD(status->statusPtr,
			    &hd, &(phmdVS->phmd[j*length+k]) ), status);
  }

  TRY( LALHOUGHIntegrHD2HT(status->statusPtr, ht, &hd), status);

  /* Free memory and exit */
  LALFree(hd.map);

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}



/**
 * This function computes the corresponding frequency bin of
 * a  \c phmd <tt>UINT8 *fBinMap</tt> for a given  intrinsic
 * frequency bin of a source <tt>UINT8 *f0Bin</tt>, and information regarding the
 * time and the residual spin down parameters HOUGHResidualSpinPar *rs.
 */
void LALHOUGHComputeFBinMap (LALStatus             *status,
			     UINT8                 *fBinMap,
			     UINT8                 *f0Bin,
			     HOUGHResidualSpinPar  *rs)
{

  UINT4    i;
  INT4    shiftFBin;
  REAL8   shiftF;

  UINT4   spinOrder;
  REAL8   *spinF;
  REAL8   timeDiff;    /*  T(t)-T(t0) */
  REAL8   timeDiffProd;
  REAL8   deltaF;  /*  df=1/TCOH  */
  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (fBinMap, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (f0Bin,   status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (rs,      status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  /*   Make sure the Input/Output pointers are not the same */
  ASSERT (fBinMap != f0Bin, status, LALHOUGHH_ESAME, LALHOUGHH_MSGESAME);

  shiftFBin = 0;
  shiftF = 0.0;

  spinOrder = rs->spinRes.length;

  if(spinOrder){
    ASSERT (rs->spinRes.data , status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
    timeDiff = rs->timeDiff;
    timeDiffProd =  timeDiff;

    deltaF = rs->deltaF;
    spinF = rs->spinRes.data;

    for (i=0; i<spinOrder; ++i ){
      shiftF += spinF[i] * timeDiffProd;
      timeDiffProd *= timeDiff;
    }
    shiftFBin = rint( shiftF/deltaF) ; /* positive or negative */
  }

  *fBinMap = *f0Bin + shiftFBin;

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}





/**
 * Calculates the total hough map for a given trajectory in the
 * time-frequency plane and a set of partial hough map derivatives allowing
 * each PHMD to have a different weight factor to account for varying
 * sensitivity at different sky-locations.
 */

void LALHOUGHConstructHMT_W (LALStatus                  *status,	/**< pointer to LALStatus structure */
			     HOUGHMapTotal              *ht, 		/**< The output hough map */
			     UINT8FrequencyIndexVector  *freqInd, 	/**< time-frequency trajectory */
			     PHMDVectorSequence         *phmdVS 	/**< set of partial hough map derivatives */)
{


  UINT4    k,j;
  UINT4    breakLine;
  UINT4    nfSize;    /* number of different frequencies */
  UINT4    length;    /* number of elements for each frequency */
  UINT8    fBinMin;   /* present minimum frequency bin */
  INT8     fBin;      /* present frequency bin */
  UINT2    xSide,ySide;

  HOUGHMapDeriv hd; /* the Hough map derivative */

  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (phmdVS,  status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (ht,      status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (freqInd, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  ASSERT (phmdVS->phmd,  status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (freqInd->data, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  /* Make sure there is no size mismatch */
  ASSERT (freqInd->length == phmdVS->length, status,
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  ASSERT (freqInd->deltaF == phmdVS->deltaF, status,
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  /* -------------------------------------------   */

  /* Make sure there are elements  */
  ASSERT (phmdVS->length, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  ASSERT (phmdVS->nfSize, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  /* -------------------------------------------   */

   /* Make sure the ht map contains some pixels */
  ASSERT (ht->xSide, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  ASSERT (ht->ySide, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);

  length = phmdVS->length;
  nfSize = phmdVS->nfSize;

  fBinMin = phmdVS->fBinMin; /* initial frequency value  od the cilinder*/

  breakLine = phmdVS->breakLine;

  /* number of physical pixels */
  xSide = ht->xSide;
  ySide = ht->ySide;

  /* Make sure initial breakLine is in [0,nfSize)  */
  ASSERT ( breakLine < nfSize, status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);

  /* -------------------------------------------   */

  /* Initializing  hd map and memory allocation */
  hd.xSide = xSide;
  hd.ySide = ySide;
  hd.map = (HoughDT *)LALMalloc(ySide*(xSide+1)*sizeof(HoughDT));
  if (hd. map == NULL) {
    ABORT( status, LALHOUGHH_EMEM, LALHOUGHH_MSGEMEM);
  }

  /* -------------------------------------------   */

  TRY( LALHOUGHInitializeHD(status->statusPtr, &hd), status);
  for ( k=0; k<length; ++k ){
    /* read the frequency index and make sure is in the proper interval*/
    fBin =freqInd->data[k] -fBinMin;

    ASSERT ( fBin < nfSize, status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);
    ASSERT ( fBin >= 0,     status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);

    /* find index */
    j = (fBin + breakLine) % nfSize;

    /* Add the corresponding PHMD to HD */
    TRY( LALHOUGHAddPHMD2HD_W(status->statusPtr,
			    &hd, &(phmdVS->phmd[j*length+k]) ), status);
  }

  TRY( LALHOUGHIntegrHD2HT(status->statusPtr, ht, &hd), status);

  /* Free memory and exit */
  LALFree(hd.map);

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}



/**
 * Adds weight factors for set of partial hough map derivatives -- the
 * weights must be calculated outside this function.
 */

void LALHOUGHWeighSpacePHMD  (LALStatus            *status,	/**< pointer to LALStatus structure */
			      PHMDVectorSequence   *phmdVS, 	/**< partial hough map derivatives */
			      REAL8Vector *weightV 		/**< vector of weights */)
{

  UINT4    k,j;
  UINT4    nfSize;    /* number of different frequencies */
  UINT4    length;    /* number of elements for each frequency */
  /* --------------------------------------------- */

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);


  /*   Make sure the arguments are not NULL: */
  ASSERT (phmdVS, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (weightV,    status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  ASSERT (phmdVS->phmd, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (weightV->data,status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  /* Make sure there is no size mismatch */
  ASSERT (weightV->length == phmdVS->length, status,
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  /* -------------------------------------------   */

  /* Make sure there are elements to be computed*/
  ASSERT (phmdVS->length, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  ASSERT (phmdVS->nfSize, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);


  length = phmdVS->length;
  nfSize = phmdVS->nfSize;

  /* weigh the phmds according to weightV */
  for ( k=0; k<length; ++k ){
    for ( j=0; j<  nfSize; ++j ){
      phmdVS->phmd[ j*length+k ].weight = (HoughDT)weightV->data[k];
    }
  }


  DETATCHSTATUSPTR (status);
   /* normal exit */
  RETURN (status);
}



/** Initializes weight factors to unity */

void LALHOUGHInitializeWeights  (LALStatus  *status,	/**< pointer to LALStatus structure */
				REAL8Vector *weightV 	/**< vector of weights */)
{

  UINT4 j, length;

   /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);


  /*   Make sure the arguments are not NULL: */
  ASSERT (weightV, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (weightV->data,status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  /* Make sure there are elements to be computed*/
  ASSERT (weightV->length, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);

  length = weightV->length;

  for (j=0; j<length; j++) {
    weightV->data[j] = 1.0;
  }

  DETATCHSTATUSPTR (status);
   /* normal exit */
  RETURN (status);
}




/** Normalizes weight factors so that their sum is N */

void LALHOUGHNormalizeWeights  (LALStatus  *status,	/**< pointer to LALStatus structure */
				REAL8Vector *weightV 	/**< vector of weights */)
{

  UINT4 j, length;
  REAL8 sum;

   /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);


  /*   Make sure the arguments are not NULL: */
  ASSERT (weightV, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (weightV->data,status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  /* Make sure there are elements to be computed*/
  ASSERT (weightV->length, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);

  length = weightV->length;

  /* calculate sum of weights */
  sum = 0.0;
  for (j=0; j<length; j++) {
    sum += weightV->data[j];
  }

  /* normalize weights */
  for (j=0; j<length; j++) {
    weightV->data[j] *= length/sum;
  }

  DETATCHSTATUSPTR (status);
   /* normal exit */
  RETURN (status);
}



/**
 * Computes weight factors arising from amplitude modulation -- it multiplies
 * an existing weight vector
 */
void LALHOUGHComputeAMWeights  (LALStatus          *status,
				REAL8Vector        *weightV,
				LIGOTimeGPSVector  *timeV,
				LALDetector        *detector,
				EphemerisData      *edat,
				REAL8              alpha,
				REAL8              delta)
{

  UINT4 length, j;

  /* amplitude modulation stuff */
  REAL4Vector *aVec, *bVec;
  REAL8 A, B, a, b;
  AMCoeffs amc;
  AMCoeffsParams *amParams;
  EarthState earth;
  BarycenterInput baryinput;

   /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);


  /*   Make sure the arguments are not NULL: */
  ASSERT (weightV, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (timeV, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (detector, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (edat, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);

  ASSERT (weightV->data,status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (timeV->data,status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  /* Make sure there is no size mismatch */
  ASSERT (weightV->length == timeV->length, status, LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  /* -------------------------------------------   */

  /* Make sure there are elements to be computed*/
  ASSERT (timeV->length, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);

  length = timeV->length;


  /* detector location */
  baryinput.site.location[0] = detector->location[0]/LAL_C_SI;
  baryinput.site.location[1] = detector->location[1]/LAL_C_SI;
  baryinput.site.location[2] = detector->location[2]/LAL_C_SI;
  baryinput.dInv = 0.e0;

  /* alpha and delta must come from the skypatch */
  /* for now set it to something arbitrary */
  /* baryinput.alpha = 0.0; */
  /* baryinput.delta = 0.0; */

  baryinput.alpha = alpha;
  baryinput.delta = delta;

  /* Allocate space for amParams stucture */
  /* Here, amParams->das is the Detector and Source info */
  amParams = (AMCoeffsParams *)LALMalloc(sizeof(AMCoeffsParams));
  if (amParams == NULL) {
    ABORT( status, LALHOUGHH_EMEM, LALHOUGHH_MSGEMEM);
  }

  amParams->das = (LALDetAndSource *)LALMalloc(sizeof(LALDetAndSource));
  if (amParams->das == NULL) {
    ABORT( status, LALHOUGHH_EMEM, LALHOUGHH_MSGEMEM);
  }

  amParams->das->pSource = (LALSource *)LALMalloc(sizeof(LALSource));
  if (amParams->das->pSource == NULL) {
    ABORT( status, LALHOUGHH_EMEM, LALHOUGHH_MSGEMEM);
  }

  /* Fill up AMCoeffsParams structure */
  amParams->baryinput = &baryinput;
  amParams->earth = &earth;
  amParams->edat = edat;
  amParams->das->pDetector = detector;
  /* make sure alpha and delta are correct */
  amParams->das->pSource->equatorialCoords.latitude = baryinput.delta;
  amParams->das->pSource->equatorialCoords.longitude = baryinput.alpha;
  amParams->das->pSource->orientation = 0.0;
  amParams->das->pSource->equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  amParams->polAngle = amParams->das->pSource->orientation ; /* These two have to be the same!!*/

  /* allocate memory for a[i] and b[i] */
  /*   TRY( LALSCreateVector( status, &(amc.a), length), status); */
  /*   TRY( LALSCreateVector( status, &(amc.b), length), status); */

  amc.a = NULL;
  amc.a = (REAL4Vector *)LALMalloc(sizeof(REAL4Vector));
  if (amc.a == NULL) {
    ABORT( status, LALHOUGHH_EMEM, LALHOUGHH_MSGEMEM);
  }

  amc.a->length = length;
  amc.a->data = (REAL4 *)LALMalloc(length*sizeof(REAL4));
  if (amc.a->data == NULL) {
    ABORT( status, LALHOUGHH_EMEM, LALHOUGHH_MSGEMEM);
  }

  amc.b = NULL;
  amc.b = (REAL4Vector *)LALMalloc(sizeof(REAL4Vector));
  if (amc.b == NULL) {
    ABORT( status, LALHOUGHH_EMEM, LALHOUGHH_MSGEMEM);
  }

  amc.b->length = length;
  amc.b->data = (REAL4 *)LALMalloc(length*sizeof(REAL4));
  if (amc.b->data == NULL) {
    ABORT( status, LALHOUGHH_EMEM, LALHOUGHH_MSGEMEM);
  }

  TRY (LALComputeAM( status->statusPtr, &amc, timeV->data, amParams), status);
  aVec = amc.a; /* a and b are as defined in JKS */
  bVec = amc.b;
  A = amc.A; /* note A is twice average of a[i]^2 */
  B = amc.B; /* note B is twice average of b[i]^2 */

  for(j=0; j<length; j++){
    a = aVec->data[j];
    b = bVec->data[j];
    weightV->data[j] *= 2.0*(a*a + b*b)/(A+B);
  }

  /* normalize weights */
  TRY( LALHOUGHNormalizeWeights( status->statusPtr, weightV ), status);

  /*   TRY( LALSDestroyVector( status, &(amc.a)), status); */
  /*   TRY( LALSDestroyVector( status, &(amc.b)), status); */
  LALFree(amc.a->data);
  LALFree(amc.b->data);
  LALFree(amc.a);
  LALFree(amc.b);

  LALFree(amParams->das->pSource);
  LALFree(amParams->das);
  LALFree(amParams);

  DETATCHSTATUSPTR (status);
   /* normal exit */
  RETURN (status);
}




/**
 * Computes weight factors arising from amplitude modulation -- it multiplies
 * an existing weight vector
 */
void LALHOUGHComputeMultiIFOAMWeights  (LALStatus          *status,
					REAL8Vector        *weightV,
					SFTCatalog         *catalog,
					EphemerisData      *edat,
					REAL8              UNUSED alpha,
					REAL8              UNUSED delta)
{

  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (weightV, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (catalog, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (edat, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);

  ASSERT (weightV->data,status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (catalog->data,status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);

  /* Make sure there is no size mismatch */
  ASSERT (weightV->length == catalog->length, status, LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);

  (void)weightV;
  (void)catalog;
  (void)edat;

  DETATCHSTATUSPTR (status);
   /* normal exit */
  RETURN (status);
}



/*@}*/
