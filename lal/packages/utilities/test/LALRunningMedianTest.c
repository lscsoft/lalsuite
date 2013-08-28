/*
*  Copyright (C) 2007 Bernd Machenschalk, David Chin, Jolien Creighton
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

#include <math.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALMalloc.h>
#include <lal/SeqFactories.h>
#include <lal/PrintVector.h>
#include <lal/LALRunningMedian.h>


/**
 * \file
 * \ingroup LALRunningMedian_h
 * \author B. Machenschalk
 *
 * \brief Program to test running median function
 *
 * \heading{Usage}
 * \code
 * LALRunningMedianTest [length blocksize [lalDebugLevel]]
 * \endcode
 *
 * \heading{Description}
 *
 * This program test the LALRunningMedian functions
 * First the proper function of the input checks is tested.
 * Then it reads an array size and a block size from the command
 * line, fills an array of the given size with random numbers,
 * computes medians of all blocks with blocksize using the
 * LALRunningMedian functions and compares the results against
 * inividually calculated medians. The test is repeated with
 * blocksize - 1 (to check for even/odd errors).
 * The default values for array length and window
 * width are 1024 and 512.
 * If a value for lalDebugLevel is given, the program
 * outputs the values of the input and median arrays
 * to files, using the PrintVector functions
 * from the support package.
 *
 */

/**\name Error Codes */
/*@{*/
#define LALRUNNINGMEDIANTESTC_ENOM 0		/**< Nominal exit */
#define LALRUNNINGMEDIANTESTC_EARG 1		/**< Error parsing command-line arguments */
#define LALRUNNINGMEDIANTESTC_ESUB 2		/**< Subroutine returned error */
#define LALRUNNINGMEDIANTESTC_EALOC 3		/**< Could not allocate data space */
#define LALRUNNINGMEDIANTESTC_EFALSE 4		/**< Medians mismatch */
#define LALRUNNINGMEDIANTESTC_EERR 5		/**< Subroutine returned wrong or no error */
/*@}*/

/** \cond DONT_DOXYGEN */
#define LALRUNNINGMEDIANTESTC_MSGENOM "Nominal exit"
#define LALRUNNINGMEDIANTESTC_MSGEARG "Error parsing command-line arguments"
#define LALRUNNINGMEDIANTESTC_MSGESUB "Subroutine returned error"
#define LALRUNNINGMEDIANTESTC_MSGEALOC "Could not allocate data space"
#define LALRUNNINGMEDIANTESTC_MSGEFALSE "Medians mismatch"
#define LALRUNNINGMEDIANTESTC_MSGEERR "Subroutine returned wrong or no error"


/* Declare lalDebugLevel */

/* global program name */
char*argv0;

/* A local macro for printing error messages */
#define EXIT( code, program, message )                                \
  if ( 1 ) {                                                          \
    if (( lalDebugLevel & LALERROR ) && (code))                       \
      LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"\
                     "        %s\n", (code), (program), __FILE__,     \
                     __LINE__, "$Id$", (message) );    \
    else if ( lalDebugLevel & LALINFO )                               \
      LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"    \
                     "        %s\n", (program), __FILE__, __LINE__,   \
                     "$Id$", (message) );              \
    return (code);                                                    \
  } else (void)(0)

#if 0
#define TOLERANCE ( 10.0 * LAL_REAL8_EPS )
#define compare_double( x, y ) \
( ( y == 0 ? ( x == 0 ? 0 : fabs((x-y)/x) ) : fabs((x-y)/y) ) > TOLERANCE )
#endif


/* prototypes */

int compare_double( double x, double y );
int compare_single( float x, float y );
static int rngmed_sortindex(const void *elem1, const void *elem2);
int testDRunningMedian(LALStatus *stat, REAL8Sequence *input, UINT4 length,
		       LALRunningMedianPar param, BOOLEAN verbose, BOOLEAN bmimpl);
int testSRunningMedian(LALStatus *stat, REAL4Sequence *input, UINT4 length,
		       LALRunningMedianPar param, BOOLEAN verbose, BOOLEAN bmimpl);


struct rngmed_val_index {
  REAL8 data;
  UINT8 index;
};


int compare_double( double x, double y )
{
  double diff;
  double denom;
  if ( x < LAL_REAL8_MIN )
  {
    if ( y < LAL_REAL8_MIN )
      return 0;
    else
      denom = fabs( y );
  }
  else
  {
    denom = fabs( x );
  }
  diff = fabs( x - y );
  if ( diff / denom > 10.0 * LAL_REAL8_EPS )
  {
    return 1;
  }
  return 0;
}


int compare_single( float x, float y )
{
  float diff;
  float denom;
  if ( x < LAL_REAL4_MIN )
  {
    if ( y < LAL_REAL4_MIN )
      return 0;
    else
      denom = fabs( y );
  }
  else
  {
    denom = fabs( x );
  }
  diff = fabs( x - y );
  if ( diff / denom > 10.0 * LAL_REAL4_EPS )
  {
    return 1;
  }
  return 0;
}

static int rngmed_sortindex(const void *elem1, const void *elem2){
  /*Used in running qsort*/

  const struct rngmed_val_index *A = elem1;
  const struct rngmed_val_index *B = elem2;
  double data1, data2;

  data1=A->data;
  data2=B->data;
  if (data1 < data2)
    return -1;
  else if (data1==data2)
    return 0;
  else
    return 1;
}


int testDRunningMedian(LALStatus *stat, REAL8Sequence *input, UINT4 length,
		       LALRunningMedianPar param, BOOLEAN verbose, BOOLEAN bmimpl) {
/* Test the LALDRunningMedian (REAL8Sequence) function by
   comparing the reults to individually calculated medians */

  REAL8 median;
  REAL8Sequence *medians=NULL;
  struct rngmed_val_index *index_block;
  UINT4 i,k;

  /* create medians vector */
  LALDCreateVector( stat, &medians, input->length - param.blocksize + 1 );
  if( stat->statusCode){
    printf("ERROR: LALDCreateVector returned status %d\n",stat->statusCode);
    EXIT( LALRUNNINGMEDIANTESTC_ESUB, argv0, LALRUNNINGMEDIANTESTC_MSGESUB );
  }

  /* call running median */
  if (bmimpl)
    LALDRunningMedian2( stat, medians, input, param );
  else
    LALDRunningMedian( stat, medians, input, param );
  if ( stat->statusCode ) {
    printf("ERROR: LALDRunningMedian returned status %d\n",stat->statusCode);
    EXIT( LALRUNNINGMEDIANTESTC_ESUB, argv0, LALRUNNINGMEDIANTESTC_MSGESUB );
  }

  /* write the vectors if verbose */
  if ( verbose ) {
    LALDPrintVector(input);
    LALDPrintVector(input);
  }

  /* allocate memory for calculating single medians */
  index_block = (struct rngmed_val_index *)LALCalloc(param.blocksize, sizeof(struct rngmed_val_index));
  if(!index_block) {
      EXIT( LALRUNNINGMEDIANTESTC_EALOC, argv0, LALRUNNINGMEDIANTESTC_MSGEALOC );
  }

  /* compare all medians */
  for(i=0;i<length-param.blocksize+1;i++) {

    /* prepare array for sort */
    for(k=0;k<param.blocksize;k++){
      index_block[k].data=input->data[k+i];
      index_block[k].index=k;
    }

    /* sort */
    qsort(index_block, param.blocksize, sizeof(struct rngmed_val_index),rngmed_sortindex);

    /* find median */
    if(param.blocksize%2==1)
      median = index_block[(param.blocksize-1)/2].data;
    else
      median = (index_block[param.blocksize/2-1].data+index_block[param.blocksize/2].data)/2;

    /* compare results */
    if(compare_double(median,medians->data[i])) {
      printf("ERROR: index:%d median:% 22.15e running median:% 22.15e mismatch:% 22.15e\n",
             i, median, medians->data[i], median - medians->data[i]);
      LALFree(index_block);
      LALFree(medians);
      LALFree(input);
      EXIT( LALRUNNINGMEDIANTESTC_EFALSE, argv0, LALRUNNINGMEDIANTESTC_MSGEFALSE );
    }
  }

  LALFree(index_block);
  LALDDestroyVector(stat,&medians);
  if ( stat->statusCode ) {
    printf("ERROR: LALDestroyVector returned status %d\n",stat->statusCode);
    EXIT( LALRUNNINGMEDIANTESTC_ESUB, argv0, LALRUNNINGMEDIANTESTC_MSGESUB );
  }
  return(0);
}




int testSRunningMedian(LALStatus *stat, REAL4Sequence *input, UINT4 length,
		       LALRunningMedianPar param, BOOLEAN verbose, BOOLEAN bmimpl) {
/* Test the LALSRunningMedian (REAL4Sequence) function by
   comparing the reults to individually calculated medians */

  REAL4 median;
  REAL4Sequence *medians=NULL;
  struct rngmed_val_index *index_block;
  UINT4 i,k;

  /* create medians vector */
  LALSCreateVector( stat, &medians,  input->length - param.blocksize + 1 );
  if( stat->statusCode){
    printf("ERROR: LALDCreateVector returned status %d\n",stat->statusCode);
    EXIT( LALRUNNINGMEDIANTESTC_ESUB, argv0, LALRUNNINGMEDIANTESTC_MSGESUB );
  }

  /* call running median */
  if (bmimpl)
    LALSRunningMedian2( stat, medians, input, param );
  else
    LALSRunningMedian( stat, medians, input, param );
  if ( stat->statusCode ) {
    printf("ERROR: LALRunningMedian returned status %d\n",stat->statusCode);
    EXIT( LALRUNNINGMEDIANTESTC_ESUB, argv0, LALRUNNINGMEDIANTESTC_MSGESUB );
  }

  /* write the vectors if verbose */
  if ( verbose ) {
    LALSPrintVector(input);
    LALSPrintVector(input);
  }

  /* allocate memory for calculating single medians */
  index_block = (struct rngmed_val_index *)LALCalloc(param.blocksize, sizeof(struct rngmed_val_index));
  if(!index_block) {
      EXIT( LALRUNNINGMEDIANTESTC_EALOC, argv0, LALRUNNINGMEDIANTESTC_MSGEALOC );
  }

  /* compare all medians */
  for(i=0;i<length-param.blocksize+1;i++) {

    /* prepare array for sort */
    for(k=0;k<param.blocksize;k++){
      index_block[k].data=input->data[k+i];
      index_block[k].index=k;
    }

    /* sort */
    qsort(index_block, param.blocksize, sizeof(struct rngmed_val_index),rngmed_sortindex);

    /* find median */
    if(param.blocksize%2==1)
      median = index_block[(param.blocksize-1)/2].data;
    else
      median = (index_block[param.blocksize/2-1].data+index_block[param.blocksize/2].data)/2;

    /* compare results */
    if(compare_single(median,medians->data[i])) {
      printf("ERROR: index:%d median:%f running median:%f mismatch\n", i, median, medians->data[i]);
      LALFree(index_block);
      LALFree(medians);
      LALFree(input);
      EXIT( LALRUNNINGMEDIANTESTC_EFALSE, argv0, LALRUNNINGMEDIANTESTC_MSGEFALSE );
    }
  }

  LALFree(index_block);
  LALSDestroyVector(stat,&medians);
  if ( stat->statusCode ) {
    printf("ERROR: LALDestroyVector returned status %d\n",stat->statusCode);
    EXIT( LALRUNNINGMEDIANTESTC_ESUB, argv0, LALRUNNINGMEDIANTESTC_MSGESUB );
  }
  return(0);
}





/**************
 **** MAIN ****
 **************/


int main( int argc, char **argv )
{
  LALStatus stat;
  UINT4 blocksize = 512, length = 1024;
  REAL4Sequence *input4=NULL;
  REAL8Sequence *input8=NULL;
  LALRunningMedianPar param;
  UINT4 i;
  BOOLEAN verbose = 0;


  /* set global program name */
  argv0 = argv[0];

  /* init status pointer */
  memset(&stat, 0, sizeof(LALStatus));

  /* init param structure */
  memset(&param, 0, sizeof(param));


  /* Parse input line. */
  if ( argc >= 3 ) {
    length = atoi(argv[1]);
    blocksize = atoi(argv[2]);
  }
  if ( argc == 4 ) {
    verbose = 1;
  }
  if (blocksize <= 3){
    fprintf(stderr,"blocksize must be >3\n");
    EXIT( LALRUNNINGMEDIANTESTC_EARG, argv0, LALRUNNINGMEDIANTESTC_MSGEARG );
  }
  if (blocksize > length){
    fprintf(stderr,"blocksize must be <= length\n");
    EXIT( LALRUNNINGMEDIANTESTC_EARG, argv0, LALRUNNINGMEDIANTESTC_MSGEARG );
  }

  /* create test input */
  LALSCreateVector( &stat, &input4, length );
  if( stat.statusCode ) {
      EXIT( LALRUNNINGMEDIANTESTC_EALOC, argv0, LALRUNNINGMEDIANTESTC_MSGEALOC );
  }
  LALDCreateVector( &stat, &input8, length );
  if( stat.statusCode ) {
      EXIT( LALRUNNINGMEDIANTESTC_EALOC, argv0, LALRUNNINGMEDIANTESTC_MSGEALOC );
  }
  for(i=0;i<length;i++)
    input4->data[i] = (input8->data[i] = (double)rand()/(double)RAND_MAX);


  /* test error conditions */

#ifndef LAL_NDEBUG
  REAL4Sequence *medians4=NULL;
  REAL8Sequence *medians8=NULL;

    if ( ! lalNoDebug )
        {

  /* create median vectors */
  LALDCreateVector( &stat, &medians8, length - blocksize + 1 );
  if( stat.statusCode ) {
      EXIT( LALRUNNINGMEDIANTESTC_EALOC, argv0, LALRUNNINGMEDIANTESTC_MSGEALOC );
  }
  LALSCreateVector( &stat, &medians4, length - blocksize + 1 );
  if( stat.statusCode ) {
      EXIT( LALRUNNINGMEDIANTESTC_EALOC, argv0, LALRUNNINGMEDIANTESTC_MSGEALOC );
  }

  /* blocksize = 0 checks */

  memset(&stat, 0, sizeof(LALStatus));
  LALSRunningMedian(&stat,medians4,input4,param);
  if(stat.statusCode == LALRUNNINGMEDIANH_EZERO) {
    printf("  PASS: LALSRunningMedian blocksize =0 results in error\n");
  } else {
    EXIT( LALRUNNINGMEDIANTESTC_EERR, argv0, LALRUNNINGMEDIANTESTC_MSGEERR );
  }

  memset(&stat, 0, sizeof(LALStatus));
  LALDRunningMedian(&stat,medians8,input8,param);
  if(stat.statusCode == LALRUNNINGMEDIANH_EZERO) {
    printf("  PASS: LALDRunningMedian blocksize =0 results in error\n");
  } else {
    EXIT( LALRUNNINGMEDIANTESTC_EERR, argv0, LALRUNNINGMEDIANTESTC_MSGEERR );
  }


  /* blocksize = 2 checks */

  param.blocksize = 2;

  memset(&stat, 0, sizeof(LALStatus));
  LALSRunningMedian(&stat,medians4,input4,param);
  if(stat.statusCode == LALRUNNINGMEDIANH_EZERO) {
    printf("  PASS: LALSRunningMedian blocksize =2 results in error\n");
  } else {
    EXIT( LALRUNNINGMEDIANTESTC_EERR, argv0, LALRUNNINGMEDIANTESTC_MSGEERR );
  }

  memset(&stat, 0, sizeof(LALStatus));
  LALDRunningMedian(&stat,medians8,input8,param);
  if(stat.statusCode == LALRUNNINGMEDIANH_EZERO) {
    printf("  PASS: LALDRunningMedian blocksize =2 results in error\n");
  } else {
    EXIT( LALRUNNINGMEDIANTESTC_EERR, argv0, LALRUNNINGMEDIANTESTC_MSGEERR );
  }


  /* blocksize too large checks */

  param.blocksize = length+1;

  memset(&stat, 0, sizeof(LALStatus));
  LALSRunningMedian(&stat,medians4,input4,param);
  if(stat.statusCode == LALRUNNINGMEDIANH_ELARGE) {
    printf("  PASS: LALSRunningMedian blocksize > input length results in error\n");
  } else {
    EXIT( LALRUNNINGMEDIANTESTC_EERR, argv0, LALRUNNINGMEDIANTESTC_MSGEERR );
  }

  memset(&stat, 0, sizeof(LALStatus));
  LALDRunningMedian(&stat,medians8,input8,param);
  if(stat.statusCode == LALRUNNINGMEDIANH_ELARGE) {
    printf("  PASS: LALDRunningMedian blocksize > input length results in error\n");
  } else {
    EXIT( LALRUNNINGMEDIANTESTC_EERR, argv0, LALRUNNINGMEDIANTESTC_MSGEERR );
  }

        } /* if ( ! lalNoDebug ) */
#endif /* LAL_NDEBUG */


  /* now set the blocksize corretly for the rest of the program */
  param.blocksize = blocksize;

#ifndef LAL_NDEBUG
    if ( ! lalNoDebug )
        {

  /* NULL pointer input check */

  memset(&stat, 0, sizeof(LALStatus));
  LALSRunningMedian(&stat,medians4,NULL,param);
  if(stat.statusCode == LALRUNNINGMEDIANH_ENULL) {
    printf("  PASS: LALSRunningMedian NULL input results in error\n");
  } else {
    EXIT( LALRUNNINGMEDIANTESTC_EERR, argv0, LALRUNNINGMEDIANTESTC_MSGEERR );
  }

  memset(&stat, 0, sizeof(LALStatus));
  LALDRunningMedian(&stat,medians8,NULL,param);
  if(stat.statusCode == LALRUNNINGMEDIANH_ENULL) {
    printf("  PASS: LALDRunningMedian NULL input results in error\n");
  } else {
    EXIT( LALRUNNINGMEDIANTESTC_EERR, argv0, LALRUNNINGMEDIANTESTC_MSGEERR );
  }

  /* median array size checks */
  memset(&stat, 0, sizeof(LALStatus));
  LALSRunningMedian(&stat,NULL,input4,param);
  if(stat.statusCode == LALRUNNINGMEDIANH_EIMED) {
    printf("  PASS: LALSRunningMedian NULL medians results in error\n");
  } else {
    EXIT( LALRUNNINGMEDIANTESTC_EERR, argv0, LALRUNNINGMEDIANTESTC_MSGEERR );
  }

  memset(&stat, 0, sizeof(LALStatus));
  LALDRunningMedian(&stat,NULL,input8,param);
  if(stat.statusCode == LALRUNNINGMEDIANH_EIMED) {
    printf("  PASS: LALDRunningMedian NULL medians results in error\n");
  } else {
    EXIT( LALRUNNINGMEDIANTESTC_EERR, argv0, LALRUNNINGMEDIANTESTC_MSGEERR );
  }

  /* destroy median test vectors */
  LALDDestroyVector(&stat,&medians8);
  if( stat.statusCode ) {
      EXIT( LALRUNNINGMEDIANTESTC_EALOC, argv0, LALRUNNINGMEDIANTESTC_MSGEALOC );
  }
  LALSDestroyVector(&stat,&medians4);
  if( stat.statusCode ) {
      EXIT( LALRUNNINGMEDIANTESTC_EALOC, argv0, LALRUNNINGMEDIANTESTC_MSGEALOC );
  }

  /* test median vectors with wrong size */


  /* too small by one */

  LALSCreateVector( &stat, &medians4, length - blocksize );
  if( stat.statusCode ) {
      EXIT( LALRUNNINGMEDIANTESTC_EALOC, argv0, LALRUNNINGMEDIANTESTC_MSGEALOC );
  }

  memset(&stat, 0, sizeof(LALStatus));
  LALSRunningMedian(&stat,medians4,input4,param);
  if(stat.statusCode == LALRUNNINGMEDIANH_EIMED) {
    printf("  PASS: LALSRunningMedian with too small median array results in error\n");
  } else {
    EXIT( LALRUNNINGMEDIANTESTC_EERR, argv0, LALRUNNINGMEDIANTESTC_MSGEERR );
  }

  LALSDestroyVector(&stat,&medians4);
  if( stat.statusCode ) {
      EXIT( LALRUNNINGMEDIANTESTC_EALOC, argv0, LALRUNNINGMEDIANTESTC_MSGEALOC );
  }

  LALDCreateVector( &stat, &medians8, length - blocksize );
  if( stat.statusCode ) {
      EXIT( LALRUNNINGMEDIANTESTC_EALOC, argv0, LALRUNNINGMEDIANTESTC_MSGEALOC );
  }

  memset(&stat, 0, sizeof(LALStatus));
  LALDRunningMedian(&stat,medians8,input8,param);
  if(stat.statusCode == LALRUNNINGMEDIANH_EIMED) {
    printf("  PASS: LALDRunningMedian with too small median array results in error\n");
  } else {
    EXIT( LALRUNNINGMEDIANTESTC_EERR, argv0, LALRUNNINGMEDIANTESTC_MSGEERR );
  }

  LALDDestroyVector(&stat,&medians8);
  if( stat.statusCode ) {
      EXIT( LALRUNNINGMEDIANTESTC_EALOC, argv0, LALRUNNINGMEDIANTESTC_MSGEALOC );
  }


  /* too large by one */

  LALSCreateVector( &stat, &medians4, length - blocksize + 2);
  if( stat.statusCode ) {
      EXIT( LALRUNNINGMEDIANTESTC_EALOC, argv0, LALRUNNINGMEDIANTESTC_MSGEALOC );
  }

  memset(&stat, 0, sizeof(LALStatus));
  LALSRunningMedian(&stat,medians4,input4,param);
  if(stat.statusCode == LALRUNNINGMEDIANH_EIMED) {
    printf("  PASS: LALSRunningMedian with too large median array results in error\n");
  } else {
    EXIT( LALRUNNINGMEDIANTESTC_EERR, argv0, LALRUNNINGMEDIANTESTC_MSGEERR );
  }

  LALSDestroyVector(&stat,&medians4);
  if( stat.statusCode ) {
      EXIT( LALRUNNINGMEDIANTESTC_EALOC, argv0, LALRUNNINGMEDIANTESTC_MSGEALOC );
  }

  LALDCreateVector( &stat, &medians8, length - blocksize + 2);
  if( stat.statusCode ) {
      EXIT( LALRUNNINGMEDIANTESTC_EALOC, argv0, LALRUNNINGMEDIANTESTC_MSGEALOC );
  }

  memset(&stat, 0, sizeof(LALStatus));
  LALDRunningMedian(&stat,medians8,input8,param);
  if(stat.statusCode == LALRUNNINGMEDIANH_EIMED) {
    printf("  PASS: LALDRunningMedian with too large median array results in error\n");
  } else {
    EXIT( LALRUNNINGMEDIANTESTC_EERR, argv0, LALRUNNINGMEDIANTESTC_MSGEERR );
  }

  LALDDestroyVector(&stat,&medians8);
  if( stat.statusCode ) {
      EXIT( LALRUNNINGMEDIANTESTC_EALOC, argv0, LALRUNNINGMEDIANTESTC_MSGEALOC );
  }

        } /* if ( ! lalNoDebug ) */
#endif /* LAL_NDEBUG */

  /* finaly restore status after checking error conditions */
  memset(&stat, 0, sizeof(LALStatus));


  /* test normal operation */

  if(testDRunningMedian(&stat,input8,length,param,verbose,0)) {
    EXIT( LALRUNNINGMEDIANTESTC_EFALSE, argv0, LALRUNNINGMEDIANTESTC_MSGEFALSE );
  } else {
    printf("  PASS: LALDRunningMedian(%d,%d)\n",length,param.blocksize);
  }

  if(testSRunningMedian(&stat,input4,length,param,verbose,0)) {
    EXIT( LALRUNNINGMEDIANTESTC_EFALSE, argv0, LALRUNNINGMEDIANTESTC_MSGEFALSE );
  } else {
    printf("  PASS: LALSRunningMedian(%d,%d)\n",length,param.blocksize);
  }

  /* decrement the blocksize for the next two test to check for even/odd errors */
  param.blocksize--;

  if(testDRunningMedian(&stat,input8,length,param,verbose,0)) {
    EXIT( LALRUNNINGMEDIANTESTC_EFALSE, argv0, LALRUNNINGMEDIANTESTC_MSGEFALSE );
  } else {
    printf("  PASS: LALDRunningMedian(%d,%d)\n",length,param.blocksize);
  }

  if(testSRunningMedian(&stat,input4,length,param,verbose,0)) {
    EXIT( LALRUNNINGMEDIANTESTC_EFALSE, argv0, LALRUNNINGMEDIANTESTC_MSGEFALSE );
  } else {
    printf("  PASS: LALSRunningMedian(%d,%d)\n",length,param.blocksize);
  }

  if(testDRunningMedian(&stat,input8,length,param,verbose,1)) {
    EXIT( LALRUNNINGMEDIANTESTC_EFALSE, argv0, LALRUNNINGMEDIANTESTC_MSGEFALSE );
  } else {
    printf("  PASS: LALDRunningMedian2(%d,%d)\n",length,param.blocksize);
  }

  if(testSRunningMedian(&stat,input4,length,param,verbose,1)) {
    EXIT( LALRUNNINGMEDIANTESTC_EFALSE, argv0, LALRUNNINGMEDIANTESTC_MSGEFALSE );
  } else {
    printf("  PASS: LALSRunningMedian2(%d,%d)\n",length,param.blocksize);
  }

  /* decrement the blocksize for the next two test to check for even/odd errors */
  param.blocksize--;

  if(testDRunningMedian(&stat,input8,length,param,verbose,1)) {
    EXIT( LALRUNNINGMEDIANTESTC_EFALSE, argv0, LALRUNNINGMEDIANTESTC_MSGEFALSE );
  } else {
    printf("  PASS: LALDRunningMedian2(%d,%d)\n",length,param.blocksize);
  }

  if(testSRunningMedian(&stat,input4,length,param,verbose,1)) {
    EXIT( LALRUNNINGMEDIANTESTC_EFALSE, argv0, LALRUNNINGMEDIANTESTC_MSGEFALSE );
  } else {
    printf("  PASS: LALSRunningMedian2(%d,%d)\n",length,param.blocksize);
  }


  /* free dummy input memory */
  LALDDestroyVector(&stat,&input8);
  LALSDestroyVector(&stat,&input4);

  /* check for memory leaks */
  LALCheckMemoryLeaks();

  /* report status if wanted */
  if(lalDebugLevel)
    REPORTSTATUS(&stat);

  /* nominal exit */
  EXIT( LALRUNNINGMEDIANTESTC_ENOM, argv0, LALRUNNINGMEDIANTESTC_MSGENOM );
}
/** \endcond */
