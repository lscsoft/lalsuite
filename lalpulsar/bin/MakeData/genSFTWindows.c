/*
 *  Copyright (C) 2022 Karl Wette
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

/*
 * Author: K. Wette
 *
 * A simple script to print out SFT window functions
 */

#include "config.h"

#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Window.h>

/* GLOBAL VARIABLES */
extern REAL8 winFncRMS;
extern REAL8TimeSeries dataDouble;
extern REAL4TimeSeries dataSingle;

/* FUNCTION PROTOTYPES */
int WindowData( REAL8 r );
int WindowDataTukey2( void );
int WindowDataHann( void );

int main( void )
{

  const size_t NWINDOWS = 4;
  const size_t WINDOWLENGTH = 256 * 1800;

  // default command line arguments from MakeSFTs.c
  const REAL8 windowR = 0.001;

  // allocate memory
  char windownames[NWINDOWS][1024];
  REAL8Vector *windows[NWINDOWS];
  REAL8 winrms[NWINDOWS];
  for ( size_t i = 0; i < NWINDOWS; ++i ) {
    windows[i] = XLALCreateREAL8Vector( WINDOWLENGTH );
    XLAL_CHECK_MAIN( windows[i] != NULL, XLAL_ENOMEM );
    for ( size_t j = 0; j < WINDOWLENGTH; ++j ) {
      windows[i]->data[j] = 1.0;
    }
  }

  size_t w = 0;

  {
    snprintf( windownames[w], sizeof( windownames[w] ), "lalpulsar_MakeSFTs Matlab style Tukey window [windowR=%g]", windowR );
    dataDouble.data = windows[w];
    WindowData( windowR );
    winrms[w] = winFncRMS;
  }

  ++w;

  {
    snprintf( windownames[w], sizeof( windownames[w] ), "lalpulsar_MakeSFTs Hann window" );
    dataDouble.data = windows[w];
    WindowDataHann();
    winrms[w] = winFncRMS;
  }

  ++w;

  {
    REAL8 param = windowR;
    REAL8Window *win = XLALCreateTukeyREAL8Window( WINDOWLENGTH, param );
    XLAL_CHECK_MAIN( win != NULL, XLAL_EFUNC );
    snprintf( windownames[w], sizeof( windownames[w] ), "XLALCreateTukeyREAL8Window(param=%g)", param );
    for ( size_t j = 0; j < WINDOWLENGTH; ++j ) {
      windows[w]->data[j] *= win->data->data[j];
    }
    winrms[w] = sqrt( win->sumofsquares / win->data->length );
    XLALDestroyREAL8Window( win );
  }

  ++w;

  {
    REAL8Window *win = XLALCreateHannREAL8Window( WINDOWLENGTH );
    XLAL_CHECK_MAIN( win != NULL, XLAL_EFUNC );
    snprintf( windownames[w], sizeof( windownames[w] ), "XLALCreateHannREAL8Window()" );
    for ( size_t j = 0; j < WINDOWLENGTH; ++j ) {
      windows[w]->data[j] *= win->data->data[j];
    }
    winrms[w] = sqrt( win->sumofsquares / win->data->length );
    XLALDestroyREAL8Window( win );
  }

  XLAL_CHECK_MAIN( w + 1 == NWINDOWS, XLAL_EFAILED );

  // output windows
  for ( size_t i = 0; i < NWINDOWS; ++i ) {
    printf( "%s%c", windownames[i], i + 1 < NWINDOWS ? ',' : '\n' );
  }
  for ( size_t i = 0; i < NWINDOWS; ++i ) {
    printf( "%0.8f%c", winrms[i], i + 1 < NWINDOWS ? ',' : '\n' );
  }
  for ( size_t j = 0; j < WINDOWLENGTH; ++j ) {
    for ( size_t i = 0; i < NWINDOWS; ++i ) {
      printf( "%0.8f%c", windows[i]->data[j], i + 1 < NWINDOWS ? ',' : '\n' );
    }
  }

  // cleanup
  for ( size_t i = 0; i < NWINDOWS; ++i ) {
    XLALDestroyREAL8Vector( windows[i] );
  }
  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

}
