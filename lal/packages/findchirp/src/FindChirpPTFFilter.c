/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpPTFFilter.c
 *
 * Author: Brown, D. A. and Fazi, D.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpPTFFilterCV">
Author: Brown, D. A. and Fazi, D.
$Id$
</lalVerbatim>

<lalLaTeX>
\input{FindChirpPTFFilterCDoc}

\vfill{\footnotesize\input{FindChirpPTFFilterCV}}
#endif

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>

double rint(double x);

NRCSID (FINDCHIRPPTFFILTERC, "$Id$");


/* <lalVerbatim file="FindChirpPTFFilterCP"> */
void
LALFindChirpPTFFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    )
/* </lalVerbatim> */
{
  UINT4                 i, j, k, l, kmax;
  UINT4                 numPoints;
  UINT4                 deltaEventIndex;
  UINT4                 ignoreIndex;
  UINT4                 haveEvent   = 0;
  REAL4                 deltaT, sum;
  REAL4                 PTFA[5][5], PTFmatrix[25];
  REAL8                 deltaF;
  REAL4                 snrThresh      = 0;
  REAL4                *snr            = NULL;
  REAL4                *PTFBinverse[5];
  COMPLEX8             *PTFQtilde, *qtilde, *PTFq, *inputData;
  COMPLEX8Vector        qVec;
  

  INITSTATUS( status, "LALFindChirpPTFFilter", FINDCHIRPPTFFILTERC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( eventList, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( !*eventList, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the filter parameters are reasonable */
  ASSERT( params->deltaT > 0, status,
      FINDCHIRPH_EDTZO, FINDCHIRPH_MSGEDTZO );
  ASSERT( params->rhosqThresh >= 0, status,
      FINDCHIRPH_ERHOT, FINDCHIRPH_MSGERHOT );

  /* check that the fft plan exists */
  ASSERT( params->invPlan, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the workspace vectors exist */

  /* if a rhosqVec vector has been created, check we can store data in it */
  if ( params->rhosqVec ) 
  {
    ASSERT( params->rhosqVec->data->data, status, 
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
    ASSERT( params->rhosqVec->data, status, 
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* make sure that the input structure exists */
  ASSERT( input, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the input structure contains some input */
  ASSERT( input->fcTmplt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->segment, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the filter has been niitialized for the correct */
  /* approximant                                                    */
  if ( params->approximant != FindChirpPTF )
  {
    ABORT( status, FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  }

  /* make sure the approximant in the tmplt and segment agree */
  ASSERT( input->fcTmplt->tmplt.approximant == FindChirpPTF, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  ASSERT( input->segment->approximant == FindChirpPTF, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );

  /*
   *
   * point local pointers to input and output pointers
   *
   */

  /* number of points in a segment */
  numPoints = params->PTFqVec->vectorLength;

  /* workspace vectors */
  snr = params->PTFsnrVec->data;
  qtilde = params->qtildeVec->data;
  PTFq   = params->PTFqVec->data;
  qVec.length = numPoints;

  /* template and data */
  inputData = input->segment->data->data->data;
  PTFQtilde = input->fcTmplt->PTFQtilde->data;
  
  for (i=0; i<5; i++) PTFBinverse[i] = input->fcTmplt->PTFBinverse->data +
                                       i * 5;

  /* number of points and frequency cutoffs */
  deltaT = params->deltaT;
  deltaF = 1.0 / ( (REAL4) params->deltaT * (REAL4) numPoints );
  kmax = input->fcTmplt->tmplt.fFinal / deltaF < numPoints/2 ? 
         input->fcTmplt->tmplt.fFinal / deltaF : numPoints/2;


  /*
   *
   * compute viable search regions in the snrsq vector
   *
   */


  if ( input->fcTmplt->tmplt.tC <= 0 )
  {
    ABORT( status, FINDCHIRPH_ECHTZ, FINDCHIRPH_MSGECHTZ );
  }

  deltaEventIndex = (UINT4) rint( (input->fcTmplt->tmplt.tC / deltaT) + 1.0 );

  /* ignore corrupted data at start and end */
  ignoreIndex = ( input->segment->invSpecTrunc / 2 ) + deltaEventIndex;

  if ( lalDebugLevel & LALINFO )
  {
    CHAR infomsg[256];

    LALSnprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg),
        "m1 = %e, m2 = %e, chi = %e, kappa = %e "
        "=> %e seconds => %d points\n"
        "invSpecTrunc = %d => ignoreIndex = %d\n", 
        input->fcTmplt->tmplt.mass1, input->fcTmplt->tmplt.mass2, 
        input->fcTmplt->tmplt.chi, input->fcTmplt->tmplt.kappa, 
        input->fcTmplt->tmplt.tC , deltaEventIndex, 
        input->segment->invSpecTrunc, ignoreIndex );
    LALInfo( status, infomsg );
  }

  /* XXX check that we are not filtering corrupted data XXX */
  /* XXX this is hardwired to 1/4 segment length        XXX */
  if ( ignoreIndex > numPoints / 4 )
  {
    ABORT( status, FINDCHIRPH_ECRUP, FINDCHIRPH_MSGECRUP );
  }
  /* XXX reset ignoreIndex to one quarter of a segment XXX */
  ignoreIndex = numPoints / 4;

  if ( lalDebugLevel & LALINFO )
  {
    CHAR infomsg[256];

    LALSnprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg), 
        "filtering from %d to %d\n",
        ignoreIndex, numPoints - ignoreIndex );
    LALInfo( status, infomsg );
  }


  /*
   *
   * compute the PTF filter statistic
   *
   */
  

  /* clear the snr output vector */
  memset( params->PTFsnrVec->data, 0, 
      params->PTFsnrVec->length * sizeof(REAL4) );

  for ( i = 0; i < params->PTFqVec->length; ++i )
  {
    
    /* compute qtilde using data and Qtilde */

    memset( params->qtildeVec->data, 0, 
            params->qtildeVec->length * sizeof(COMPLEX8) );

    /* qtilde positive frequency, not DC or nyquist */
    for ( k = 1; k < kmax; ++k )
    {
      REAL4 r = inputData[k].re;
      REAL4 s = inputData[k].im;
      REAL4 x = PTFQtilde[i * (numPoints / 2 + 1) + k].re;
      REAL4 y = 0 - PTFQtilde[i * (numPoints / 2 + 1) + k].im; /* cplx conj */

      qtilde[k].re = r*x - s*y;
      qtilde[k].im = r*y + s*x;
    }

    qVec.data = params->PTFqVec->data + (i * params->PTFqVec->vectorLength);
    
    /* inverse fft to get q */
    LALCOMPLEX8VectorFFT( status->statusPtr, &qVec, params->qtildeVec, 
        params->invPlan );
    CHECKSTATUSPTR( status );
  }

  /* now we have PTFqVec which contains <s|Q^I_0> + i <s|Q^I_\pi/2> */

  /* construct A */
  
  for (k=0; k<numPoints; i++)
  {  
    for (i=0; i<5; i++)
    {  
      for (j=0; j<i+1; j++)
      {  
        PTFA[i][j] = PTFq[i * numPoints + k].re * PTFq[j * numPoints + k].re +
                     PTFq[i * numPoints + k].im * PTFq[j * numPoints + k].im;
      }  
    }  
 /* multiply by PTFBinverse to obtain AB^(-1) */
    sum = 0.0;
    for (i=0; i<5; i++)
    {  
      for (j=0; j<i+1; j++)
      {  
       for (l=0; l<5; l++) sum = sum + PTFA[i][l] * PTFBinverse[l][j];
       PTFmatrix[i + 5 * j] = sum;
       sum = 0.0;
      }
    }
   
    /* find max eigenvalue and store it in snr vector */
  }
 

  fprintf( stderr, "Ptf filtering data segment\n" );

  /*
   *
   * look for and cluster events in the snr vector
   *
   */

  
  /* look for an event in the filter output */
  for ( j = ignoreIndex; j < numPoints - ignoreIndex; ++j )
  {
    /* if snrsq exceeds threshold at any point */
    if ( snr[j] > snrThresh )
    {
      haveEvent = 1;        /* mark segment to have events    */
      break;
    }
  }

  /* search the SNR vector for events */
  /* process events in the filter output */
  if ( haveEvent )
  {
#if 0
    LALFindChirpClusterEvents( status->statusPtr, eventList, input,
        params, q, kmax, numPoints, ignoreIndex, 
        norm, modqsqThresh, chisqThreshFac, numChisqBins, searchName );
    CHECKSTATUSPTR( status );
#endif
    fprintf( stderr, "found events!\n" );
  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
