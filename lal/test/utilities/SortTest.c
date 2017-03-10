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
 * ### Usage ###
 *
 * \code
 * SortTest
 * \endcode
 *
 * ### Description ###
 *
 * This test program creates rank and index arrays for an unordered list
 * of numbers, and then sorts the list.  The data for the list are
 * generated randomly, and the output is to \c stdout if <tt>-v</tt> is
 * specified (unless redirected).  \c SortTest returns 0 if it executes
 * successfully, and 1 if any of the subroutines fail.
 *
 * ### Exit codes ###
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
#include <lal/Sort.h>


static void makedata( int nobj, int **data, int **sort, int **indx, int **rank )
{
  int i;

  *data = malloc(nobj*sizeof(**data));
  *sort = malloc(nobj*sizeof(**data));
  *indx = malloc(nobj*sizeof(**indx));
  *rank = malloc(nobj*sizeof(**rank));

  for ( i = 0; i < nobj; ++i )
    (*sort)[i] = (*data)[i] = rand() % 100;
}


static void freedata( int *data, int *sort, int *indx, int *rank)
{
  free( data );
  free( sort );
  free( indx );
  free( rank );
}


static int compar( void *p, const void *a, const void *b )
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


static int check( int *data, int *sort, int nobj, int ascend )
{
  int i, j;
  /* make sure result is in correct order */
  for ( i = 1; i < nobj; ++i )
    if ( ascend )
    {
      if ( sort[i] < sort[i-1] )
        abort();
    }
    else
    {
      if ( sort[i] > sort[i-1] )
        abort();
    }
  /* make sure every element from the original data can be found somewhere
   * in the result, and every element in the result can be found in the
   * original data */
  for ( i = 0; i < nobj; ++i )
  {
    for ( j = 0; j < nobj && sort[j] != data[i]; ++j );
    if ( j == nobj )
      abort();
  }
  for ( i = 0; i < nobj; ++i )
  {
    for ( j = 0; j < nobj && data[j] != sort[i]; ++j );
    if ( j == nobj )
      abort();
  }
  return 0;
}


int main(int argc, char **argv)
{
  int testnum;

  if (argc!=1)
      LALPrintError("%s: Incorrect arguments\n",argv[0]);

  srand(time(NULL));

  for ( testnum = 0; testnum < 200; testnum++ )
  {
    int nobj = rand() % 24;
    int ascend = rand() & 1;
    int *data;
    int *sort;
    int *indx;
    int *rank;
    int i;

    makedata( nobj, &data, &sort, &indx, &rank );

    if ( XLALHeapIndex( indx, data, nobj, sizeof(*data), &ascend, compar ) < 0 )
      abort();
    if ( XLALHeapRank( rank, data, nobj, sizeof(*data), &ascend, compar ) < 0 )
      abort();
    if ( XLALHeapSort( sort, nobj, sizeof(*data), &ascend, compar ) < 0 )
      abort();

    for ( i = 0; i < nobj; ++i )
      if ( sort[i] != data[indx[i]] )
        abort();

    check( data, sort, nobj, ascend );

    freedata( data, sort, indx, rank );
  }

  for ( testnum = 0; testnum < 200; testnum++ )
  {
    int nobj = rand() % 24;
    int ascend = rand() & 1;
    int *data;
    int *sort;
    int *indx;	/* unused for these tests */
    int *rank;	/* unused for these tests */

    makedata( nobj, &data, &sort, &indx, &rank );

    if ( XLALInsertionSort( sort, nobj, sizeof(*data), &ascend, compar ) < 0 )
      abort();

    check( data, sort, nobj, ascend );

    freedata( data, sort, indx, rank );
  }


  return 0;
}

/** \endcond */
