#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>

#include "coh_PTF.h"

void calculate_bmatrix(
  struct coh_PTF_params   *params,
  gsl_matrix *eigenvecs,
  gsl_vector *eigenvals,
  REAL4 a[LAL_NUM_IFO],
  REAL4 b[LAL_NUM_IFO],
  REAL8Array              *PTFM[LAL_NUM_IFO+1],
  UINT4 vecLength,
  UINT4 vecLengthTwo,
  UINT4 PTFMlen
  )
{
  // This function calculates the eigenvectors and eigenvalues of the
  // coherent "B" matrix. This is the matrix that appears as B^{-1} in the
  // original definition of SNR.
  UINT4 i,j,k;
  UINT4 vecLengthSquare = vecLength*vecLength;
  REAL4 zh[vecLengthSquare],sh[vecLengthSquare],yu[vecLengthSquare];
  gsl_matrix *B2 = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_eigen_symmv_workspace *matTemp = gsl_eigen_symmv_alloc (vecLengthTwo);
  /* Create and invert the Bmatrix */
  for (i = 0; i < vecLength; i++ )
  {
    for (j = 0; j < vecLength; j++ )
    {
      zh[i*vecLength+j] = 0;
      sh[i*vecLength+j] = 0;
      yu[i*vecLength+j] = 0;
      for( k = 0; k < LAL_NUM_IFO; k++)
      {
        if ( params->haveTrig[k] )
        {
          zh[i*vecLength+j] += a[k]*a[k] * PTFM[k]->data[i*PTFMlen+j];
          sh[i*vecLength+j] += b[k]*b[k] * PTFM[k]->data[i*PTFMlen+j];
          yu[i*vecLength+j] += a[k]*b[k] * PTFM[k]->data[i*PTFMlen+j];
        }
      }
    }
  }

  for (i = 0; i < vecLengthTwo; i++ )
  {
    for (j = 0; j < vecLengthTwo; j++ )
    {
      if ( i < vecLength && j < vecLength )
      {
        gsl_matrix_set(B2,i,j,zh[i*vecLength+j]);
      }
      else if ( i > (vecLength-1) && j > (vecLength-1))
      {
        gsl_matrix_set(B2,i,j,sh[(i-vecLength)*vecLength + (j-vecLength)]);
      }
      else if ( i < vecLength && j > (vecLength-1))
      {
        gsl_matrix_set(B2,i,j, yu[i*vecLength + (j-vecLength)]);
      }
      else if ( i > (vecLength-1) && j < vecLength)
      {
        gsl_matrix_set(B2,i,j,yu[j*vecLength + (i-vecLength)]);
      }
      else
        fprintf(stderr,"BUGGER! Something went wrong.");
      /*fprintf(stdout,"%f ",gsl_matrix_get(B2,i,j));*/
    }
    /*fprintf(stdout,"\n");*/
  }

  /* Here we compute the eigenvalues and eigenvectors of B2 */
  gsl_eigen_symmv (B2,eigenvals,eigenvecs,matTemp);
  gsl_eigen_symmv_free(matTemp);
  gsl_matrix_free(B2);
}


void calculate_rotated_vectors(
    struct coh_PTF_params   *params,
    COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1],
    REAL4 *u1,
    REAL4 *u2,
    REAL4 a[LAL_NUM_IFO],
    REAL4 b[LAL_NUM_IFO],
    INT4  timeOffsetPoints[LAL_NUM_IFO],
    gsl_matrix *eigenvecs,
    gsl_vector *eigenvals,
    UINT4 numPoints,
    UINT4 position,
    UINT4 vecLength,
    UINT4 vecLengthTwo)
{
  // This function calculates the coherent time series and rotates them into
  // the basis where the B matrix is the identity.
  
  UINT4 j,k;  
  REAL4 v1[vecLengthTwo],v2[vecLengthTwo];

  for ( j = 0; j < vecLengthTwo ; j++ ) /* Construct the vi vectors */
  {
    v1[j] = 0.;
    v2[j] = 0.;
    for( k = 0; k < LAL_NUM_IFO; k++)
    {
      if ( params->haveTrig[k] )
      {
        if (j < vecLength)
        {
          v1[j] += a[k] * PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]].re;
          v2[j] += a[k] * PTFqVec[k]->data[j*numPoints+position+timeOffsetPoints[k]].im;
        }
        else
        {
          v1[j] += b[k] * PTFqVec[k]->data[(j-vecLength)*numPoints+position+timeOffsetPoints[k]].re;
          v2[j] += b[k] * PTFqVec[k]->data[(j-vecLength)*numPoints+position+timeOffsetPoints[k]].im;
        }
      }
    }
  }

  /* Now we rotate the v1 and v2 to be in orthogonal basis */
  /* We can use gsl multiplication stuff to do this */
  /* BLAS stuff is stupid so we'll do it explicitly! */

  for ( j = 0 ; j < vecLengthTwo ; j++ )
  {
    u1[j] = 0.;
    u2[j] = 0.;
    for ( k = 0 ; k < vecLengthTwo ; k++ )
    {
      u1[j] += gsl_matrix_get(eigenvecs,k,j)*v1[k];
      u2[j] += gsl_matrix_get(eigenvecs,k,j)*v2[k];
    }
    u1[j] = u1[j] / (pow(gsl_vector_get(eigenvals,j),0.5));
    u2[j] = u2[j] / (pow(gsl_vector_get(eigenvals,j),0.5));
  }

}
