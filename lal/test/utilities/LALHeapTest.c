/*
 *  Copyright (C) 2016 Karl Wette
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

#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <lal/LALHeap.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

static void *new_int( int x )
{
  return memcpy( XLALMalloc( sizeof( x ) ), &x, sizeof( x ) );
}

static int cmp_ptr_int( const void *x, const void *y )
{
  return *( ( const int * ) x ) - *( ( const int * ) y );
}

static int cmp_ptr_ptr_int( const void *x, const void *y )
{
  return cmp_ptr_int( *( ( const void *const * ) x ), *( ( const void *const * ) y ) );
}

static int print_ptr_int( UNUSED void *param, const void *x )
{
  printf( " %i", *( ( const int * ) x ) );
  return XLAL_SUCCESS;
}

static int check_ptr_int( void *param, const void *x )
{
  int ***y = ( int ** * ) param;
  return *( ( const int * ) x ) == **( ( *y )++ ) ? XLAL_SUCCESS : XLAL_FAILURE;
}

static int reverse_ptr_int( void *param, void *x )
{
  *( ( int * ) x ) = *( ( int * ) param ) - *( ( int * ) x );
  return XLAL_SUCCESS;
}

int main( void )
{

  /* Create heaps for storing integers */
  LALHeap *minh = XLALHeapCreate( XLALFree, 0, -1, cmp_ptr_int );
  XLAL_CHECK_MAIN( minh != NULL, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALHeapMaxSize( minh ) == 0, XLAL_EFAILED );
  LALHeap *maxh = XLALHeapCreate( XLALFree, 0, +1, cmp_ptr_int );
  XLAL_CHECK_MAIN( maxh != NULL, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALHeapMaxSize( minh ) == 0, XLAL_EFAILED );
  LALHeap *min10h = XLALHeapCreate( XLALFree, 10, -1, cmp_ptr_int );
  XLAL_CHECK_MAIN( min10h != NULL, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALHeapMaxSize( min10h ) == 10, XLAL_EFAILED );
  LALHeap *max10h = XLALHeapCreate( XLALFree, 10, +1, cmp_ptr_int );
  XLAL_CHECK_MAIN( max10h != NULL, XLAL_EFUNC );
  XLAL_CHECK_MAIN( XLALHeapMaxSize( min10h ) == 10, XLAL_EFAILED );

  /* Check some properties of empty heaps */
  xlalErrno = 0;
  XLAL_CHECK_MAIN( XLALHeapSize( minh ) == 0 && xlalErrno == 0, XLAL_EFAILED );
  xlalErrno = 0;
  XLAL_CHECK_MAIN( XLALHeapRoot( minh ) == NULL && xlalErrno == 0, XLAL_EFAILED );

  /* Create an array of 100 random integers for input */
  int input[100];
  {
    gsl_rng *r = gsl_rng_alloc( gsl_rng_mt19937 );
    XLAL_CHECK_MAIN( r != NULL, XLAL_ESYS );
    for ( size_t i = 0; i < 100; ++i ) {
      input[i] = -37 + gsl_rng_uniform_int( r, 100 );
    }
    gsl_rng_free( r );
  }

  /* Add each integer to each heap, and to a sorted reference array */
  int *ref[100];
  for ( int i = 0; i < 100; ++i ) {
    printf( "\n----- i=%i -----\n", i );

    /* Add integer to reference array, and re-sort */
    {
      ref[i] = &input[i];
      qsort( ref, i+1, sizeof( ref[0] ), cmp_ptr_ptr_int );
      printf( "ref={" );
      for ( int j = 0; j <= i; ++j ) {
        printf( " %i", *ref[j] );
      }
      printf( " }\n" );
    }

    /* Add integer to unlimited min-heap, and check properties */
    {
      void *x = NULL;
      x = new_int( input[i] );
      XLAL_CHECK_MAIN( x != NULL, XLAL_ENOMEM );
      XLAL_CHECK_MAIN( XLALHeapAdd( minh, &x ) == XLAL_SUCCESS, XLAL_EFUNC );
      printf( "minh={" );
      XLAL_CHECK_MAIN( XLALHeapVisit( minh, print_ptr_int, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
      printf( " }" );
      XLAL_CHECK_MAIN( x == NULL, XLAL_EFAILED );
      XLAL_CHECK_MAIN( XLALHeapSize( minh ) == i+1, XLAL_EFAILED );
      XLAL_CHECK_MAIN( XLALHeapRoot( minh ) != NULL, XLAL_EFAILED );
      const int r = *( ( const int * ) XLALHeapRoot( minh ) );
      const int r_ref = *ref[0];
      XLAL_CHECK_MAIN( r == r_ref, XLAL_EFAILED, "(root) %i != %i (ref)", r, r_ref );
      printf( ", root = %i\n", r );
      {
        int **ref0 = &ref[0];
        XLAL_CHECK_MAIN( XLALHeapVisit( minh, check_ptr_int, &ref0 ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
    }

    /* Add integer to unlimited max-heap, and check properties */
    {
      void *x = NULL;
      x = new_int( input[i] );
      XLAL_CHECK_MAIN( x != NULL, XLAL_ENOMEM );
      XLAL_CHECK_MAIN( XLALHeapAdd( maxh, &x ) == XLAL_SUCCESS, XLAL_EFUNC );
      printf( "maxh={" );
      XLAL_CHECK_MAIN( XLALHeapVisit( maxh, print_ptr_int, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
      printf( " }" );
      XLAL_CHECK_MAIN( x == NULL, XLAL_EFAILED );
      XLAL_CHECK_MAIN( XLALHeapSize( maxh ) == i+1, XLAL_EFAILED );
      XLAL_CHECK_MAIN( XLALHeapRoot( maxh ) != NULL, XLAL_EFAILED );
      const int r = *( ( const int * ) XLALHeapRoot( maxh ) );
      const int r_ref = *ref[i];
      XLAL_CHECK_MAIN( r == r_ref, XLAL_EFAILED, "(root) %i != %i (ref)", r, r_ref );
      printf( ", root = %i\n", r );
      {
        int **ref0 = &ref[0];
        XLAL_CHECK_MAIN( XLALHeapVisit( maxh, check_ptr_int, &ref0 ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
    }

    /* Add integer to limited min-heap, and check properties */
    {
      void *x = new_int( input[i] );
      XLAL_CHECK_MAIN( x != NULL, XLAL_ENOMEM );
      XLAL_CHECK_MAIN( XLALHeapAdd( min10h, &x ) == XLAL_SUCCESS, XLAL_EFUNC );
      printf( "min10h={" );
      XLAL_CHECK_MAIN( XLALHeapVisit( min10h, print_ptr_int, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
      printf( " }" );
      if ( i < 10 ) {
        XLAL_CHECK_MAIN( x == NULL, XLAL_EFAILED );
        XLAL_CHECK_MAIN( XLALHeapSize( min10h ) == i+1, XLAL_EFAILED );
        XLAL_CHECK_MAIN( XLALHeapRoot( min10h ) != NULL, XLAL_EFAILED );
      } else {
        XLAL_CHECK_MAIN( x != NULL, XLAL_EFAILED );
        XLALFree( x );
        XLAL_CHECK_MAIN( XLALHeapSize( min10h ) == 10, XLAL_EFAILED );
        XLAL_CHECK_MAIN( XLALHeapRoot( min10h ) != NULL, XLAL_EFAILED );
      }
      const int r = *( ( const int * ) XLALHeapRoot( min10h ) );
      const int r_ref = *ref[( i < 10 ) ? 0 : i-9];
      XLAL_CHECK_MAIN( r == r_ref, XLAL_EFAILED, "(root) %i != %i (ref)", r, r_ref );
      printf( ", root = %i\n", r );
      {
        int **ref0 = &ref[( i < 10 ) ? 0 : i-9];
        XLAL_CHECK_MAIN( XLALHeapVisit( min10h, check_ptr_int, &ref0 ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
    }

    /* Add integer to limited max-heap, and check properties */
    {
      void *x = new_int( input[i] );
      XLAL_CHECK_MAIN( x != NULL, XLAL_ENOMEM );
      XLAL_CHECK_MAIN( XLALHeapAdd( max10h, &x ) == XLAL_SUCCESS, XLAL_EFUNC );
      printf( "max10h={" );
      XLAL_CHECK_MAIN( XLALHeapVisit( max10h, print_ptr_int, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
      printf( " }" );
      if ( i < 10 ) {
        XLAL_CHECK_MAIN( x == NULL, XLAL_EFAILED );
        XLAL_CHECK_MAIN( XLALHeapSize( max10h ) == i+1, XLAL_EFAILED );
        XLAL_CHECK_MAIN( XLALHeapRoot( max10h ) != NULL, XLAL_EFAILED );
      } else {
        XLAL_CHECK_MAIN( x != NULL, XLAL_EFAILED );
        XLALFree( x );
        XLAL_CHECK_MAIN( XLALHeapSize( max10h ) == 10, XLAL_EFAILED );
        XLAL_CHECK_MAIN( XLALHeapRoot( max10h ) != NULL, XLAL_EFAILED );
      }
      const int r = *( ( const int * ) XLALHeapRoot( max10h ) );
      const int r_ref = *ref[( i < 10 ) ? i : 9];
      XLAL_CHECK_MAIN( r == r_ref, XLAL_EFAILED, "(root) %i != %i (ref)", r, r_ref );
      printf( ", root = %i\n", r );
      {
        int **ref0 = &ref[0];
        XLAL_CHECK_MAIN( XLALHeapVisit( max10h, check_ptr_int, &ref0 ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
    }

    fflush( stdout );
  }

  /* Modify unlimited heaps, and check properties */
  {
    printf( "\n----- XLALHeapModify() -----\n" );
    int param = 77;
    {
      for ( int j = 0; j < 100; ++j ) {
        *ref[j] = param - *ref[j];
      }
      qsort( ref, 100, sizeof( ref[0] ), cmp_ptr_ptr_int );
      printf( "ref={" );
      for ( int j = 0; j < 100; ++j ) {
        printf( " %i", *ref[j] );
      }
      printf( " }\n" );
    }
    {
      XLAL_CHECK( XLALHeapModify( minh, reverse_ptr_int, &param ) == XLAL_SUCCESS, XLAL_EFUNC );
      printf( "minh={" );
      XLAL_CHECK_MAIN( XLALHeapVisit( minh, print_ptr_int, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
      const int r = *( ( const int * ) XLALHeapRoot( minh ) );
      const int r_ref = *ref[0];
      XLAL_CHECK_MAIN( r == r_ref, XLAL_EFAILED, "(root) %i != %i (ref)", r, r_ref );
      printf( ", root = %i\n", r );
      {
        int **ref0 = &ref[0];
        XLAL_CHECK_MAIN( XLALHeapVisit( minh, check_ptr_int, &ref0 ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
    }
    {
      XLAL_CHECK( XLALHeapModify( maxh, reverse_ptr_int, &param ) == XLAL_SUCCESS, XLAL_EFUNC );
      printf( "maxh={" );
      XLAL_CHECK_MAIN( XLALHeapVisit( maxh, print_ptr_int, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
      const int r = *( ( const int * ) XLALHeapRoot( maxh ) );
      const int r_ref = *ref[99];
      XLAL_CHECK_MAIN( r == r_ref, XLAL_EFAILED, "(root) %i != %i (ref)", r, r_ref );
      printf( ", root = %i\n", r );
      {
        int **ref0 = &ref[0];
        XLAL_CHECK_MAIN( XLALHeapVisit( maxh, check_ptr_int, &ref0 ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
    }
    fflush( stdout );
  }

  /* Resize unlimited heaps, and check properties */
  {
    printf( "\n----- XLALHeapResize() -----\n" );
    {
      XLAL_CHECK_MAIN( XLALHeapResize( minh, 20 ) == XLAL_SUCCESS, XLAL_EFUNC );
      int **ref0 = &ref[80];
      XLAL_CHECK_MAIN( XLALHeapVisit( minh, check_ptr_int, &ref0 ) == XLAL_SUCCESS, XLAL_EFUNC );
      const int **elems = ( const int ** ) XLALHeapElements( minh );
      for ( int j = 0; j < XLALHeapSize( minh ); ++j ) {
        XLAL_CHECK_MAIN( *elems[j] == *ref[80 + j], XLAL_EFAILED, "(elems) %i != %i (ref)", *elems[j], *ref[80 + j] );
      }
      XLALFree( elems );
    }
    {
      XLAL_CHECK_MAIN( XLALHeapResize( maxh, 33 ) == XLAL_SUCCESS, XLAL_EFUNC );
      int **ref0 = &ref[0];
      XLAL_CHECK_MAIN( XLALHeapVisit( maxh, check_ptr_int, &ref0 ) == XLAL_SUCCESS, XLAL_EFUNC );
      const int **elems = ( const int ** ) XLALHeapElements( maxh );
      for ( int j = 0; j < XLALHeapSize( maxh ); ++j ) {
        XLAL_CHECK_MAIN( *elems[j] == *ref[j], XLAL_EFAILED, "(elems) %i != %i (ref)", *elems[j], *ref[j] );
      }
      XLALFree( elems );
    }
    fflush( stdout );
  }

  /* Cleanup */
  {
    printf( "\n----- cleanup -----\n" );
    XLALHeapDestroy( minh );
    XLALHeapDestroy( maxh );
    XLALHeapDestroy( min10h );
    XLALHeapDestroy( max10h );

    /* Check for memory leaks */
    LALCheckMemoryLeaks();

    fflush( stdout );
  }

  return EXIT_SUCCESS;

}
