/*
*  Copyright (C) 2007 Jolien Creighton
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
 * File Name: SortTest.c
 *
 * Author: Creighton, T. D.
 *
 *
 *-----------------------------------------------------------------------*/

/**
 * \file
 * \ingroup Sort_h
 *
 * \brief A program to test sorting routines.
 *
 * \heading{Usage}
 * \code
 * SortTest [-s seed]
 * \endcode
 *
 * \heading{Description}
 *
 * This test program creates rank and index arrays for an unordered list
 * of numbers, and then sorts the list.  The data for the list are
 * generated randomly, and the output is to \c stdout if <tt>-v</tt> is
 * specified (unless redirected).  \c SortTest returns 0 if it executes
 * successfully, and 1 if any of the subroutines fail.
 *
 * The <tt>-s</tt> option sets the seed for the random number generator; if
 * \c seed is set to zero (or if no <tt>-s</tt> option is given) then
 * the seed is taken from the processor clock.
 *
 * \heading{Exit codes}
 *
 * <table><tr><th>Code</th><th>Explanation</th></tr>
 * <tr><td>0</td><td>Success, normal exit.</td></tr>
 * <tr><td>1</td><td>Subroutine failed.</td></tr>
 * </table>
 *
 */

/** \cond DONT_DOXYGEN */

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/Sort.h>

#define NPTS 5

#define SORTTEST_ESUB 1
#define SORTTEST_MSGESUB "Subroutine returned error"

void test_xlal_routines( void );

int main(int argc, char **argv)
{
  static LALStatus stat;
  INT4         i;
  INT4         seed=0;
  INT4         verbose=0;
  INT4Vector   *lal_index=NULL;
  REAL4Vector  *data=NULL;
  RandomParams *params=NULL;

  if(argc==2&&!strcmp(argv[1],"-s"))
    seed=atoi(argv[1]);
  else if (argc!=1)
      LALPrintError("%s: Incorrect arguments\n",argv[0]);
  test_xlal_routines();

  /* Create vectors and random parameters. */
  LALI4CreateVector(&stat,&lal_index,NPTS);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return SORTTEST_ESUB;
  }
  LALSCreateVector(&stat,&data,NPTS);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return SORTTEST_ESUB;
  }
  LALCreateRandomParams(&stat,&params,seed);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return SORTTEST_ESUB;
  }

  /* Initialize data. */
  for(i=0;i<NPTS;i++){
    LALUniformDeviate(&stat,data->data+i,params);
    if(stat.statusCode){
      LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
      REPORTSTATUS(&stat);
      return SORTTEST_ESUB;
    }
    data->data[i]*=9.99;
  }

  /* Rank data; print data and ranks. */
  LALSHeapRank(&stat,lal_index,data);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return SORTTEST_ESUB;
  }
  if ( verbose )
    fprintf(stdout,  " Data    Rank\n");
  for(i=0;i<NPTS;i++)
    if ( verbose )
      fprintf(stdout," %4.2f     %2i\n",data->data[i],lal_index->data[i]);
  if ( verbose )
    fprintf(stdout,"\n");

  /* Index and sort data; print sorted data and indecies. */
  LALSHeapIndex(&stat,lal_index,data);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return SORTTEST_ESUB;
  }
  LALSHeapSort(&stat,data);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return SORTTEST_ESUB;
  }
  if ( verbose )
    fprintf(stdout,  "Sorted  Index\n");
  for(i=0;i<NPTS;i++)
    if ( verbose )
      fprintf(stdout," %4.2f     %2i\n",data->data[i],lal_index->data[i]);

  /* Free and clear. */
  LALI4DestroyVector(&stat,&lal_index);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return SORTTEST_ESUB;
  }
  LALSDestroyVector(&stat,&data);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return SORTTEST_ESUB;
  }
  LALDestroyRandomParams(&stat,&params);
  if(stat.statusCode){
    LALPrintError("%s: %s\n",argv[0],SORTTEST_MSGESUB);
    REPORTSTATUS(&stat);
    return SORTTEST_ESUB;
  }
  return 0;
}


int compar( void *p, const void *a, const void *b );
int compar( void *p, const void *a, const void *b )
{
  int x = *((const int *)a);
  int y = *((const int *)b);
  int ascend = *(int *)p;

  if ( ascend )
  {
    if ( x < y )
      return -1;
    if ( x > y )
      return 1;
    return 0;
  }

  if ( x > y )
    return -1;
  if ( x < y )
    return 1;
  return 0;
}

int check( int *data, int nobj, int ascend );
int check( int *data, int nobj, int ascend )
{
  int i;
  for ( i = 1; i < nobj; ++i )
    if ( ascend )
    {
      if ( data[i] < data[i-1] )
        abort();
    }
    else
    {
      if ( data[i] > data[i-1] )
        abort();
    }
  return 0;
}

void test_xlal_routines( void )
{
  int  nobj = 9;
  int *data;
  int *sort;
  int *indx;
  int *rank;
  int ascend;
  int code;
  int i;

  data = malloc(nobj*sizeof(*data));
  sort = malloc(nobj*sizeof(*data));
  indx = malloc(nobj*sizeof(*indx));
  rank = malloc(nobj*sizeof(*rank));

  srand(time(NULL));

  /* sort in ascending order */

  ascend = 1;

  for ( i = 0; i < nobj; ++i )
    sort[i] = data[i] = rand() % 100;

  code = XLALHeapIndex( indx, data, nobj, sizeof(*data), &ascend, compar );
  if ( code < 0 )
    abort();
  code = XLALHeapRank( rank, data, nobj, sizeof(*data), &ascend, compar );
  if ( code < 0 )
    abort();
  code = XLALHeapSort( sort, nobj, sizeof(*data), &ascend, compar );
  if ( code < 0 )
    abort();

  for ( i = 0; i < nobj; ++i )
    if ( sort[i] != data[indx[i]] )
      abort();

  check( sort, nobj, ascend );

  /* sort in descending order */

  ascend = 0;

  for ( i = 0; i < nobj; ++i )
    sort[i] = data[i] = rand() % 100;

  code = XLALHeapIndex( indx, data, nobj, sizeof(*data), &ascend, compar );
  if ( code < 0 )
    abort();
  code = XLALHeapRank( rank, data, nobj, sizeof(*data), &ascend, compar );
  if ( code < 0 )
    abort();
  code = XLALHeapSort( sort, nobj, sizeof(*data), &ascend, compar );
  if ( code < 0 )
    abort();

  for ( i = 0; i < nobj; ++i )
    if ( sort[i] != data[indx[i]] )
      abort();

  check( sort, nobj, ascend );

  free( data );
  free( sort );
  free( indx );
  free( rank );

  return;
}

/** \endcond */
