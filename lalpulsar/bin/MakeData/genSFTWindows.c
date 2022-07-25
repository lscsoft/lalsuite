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

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Window.h>

/* STRUCTURES */
struct CommandLineArgsTag {
  REAL8 HPf;              /* High pass filtering frequency */
  INT4 T;                 /* SFT duration */
  char *stringT;          /* 12/27/05 gam; string with SFT duration */
  INT4 GPSStart;
  INT4 GPSEnd;
  INT4 makeGPSDirs;        /* 12/27/05 gam; add option to make directories based on gps time */
  INT4 sftVersion;         /* 12/28/05 gam; output SFT version */
  char *commentField;      /* 12/28/05 gam; string comment for version 2 SFTs */
  BOOLEAN htdata;          /* flag that indicates we're doing h(t) data */
  BOOLEAN makeTmpFile;     /* 01/09/06 gam */
  char *FrCacheFile;       /* Frame cache file */
  char *ChannelName;
  char *IFO;               /* 01/14/07 gam */
  char *SFTpath;           /* path to SFT file location */
  char *miscDesc;          /* 12/28/05 gam; string giving misc. part of the SFT description field in the filename */
  INT4 PSSCleaning;	   /* 1=YES and 0=NO*/
  REAL8 PSSCleanHPf;       /* Cut frequency for the bilateral highpass filter. It has to be used only if PSSCleaning is YES.*/
  INT4 PSSCleanExt;        /* Extend the timeseries at the beginning before calculating the autoregressive mean */
  INT4 windowOption;       /* 12/28/05 gam; window options; 0 = no window, 1 = default = Matlab style Tukey window; 2 = make_sfts.c Tukey window; 3 = Hann window */
  REAL8 windowR;
  REAL8 overlapFraction;   /* 12/28/05 gam; overlap fraction (for use with windows; e.g., use -P 0.5 with -w 3 Hann windows; default is 1.0). */
  BOOLEAN useSingle;       /* 11/19/05 gam; use single rather than double precision */
  char *frameStructType;   /* 01/10/07 gam */
};

/* GLOBAL VARIABLES */
REAL8 winFncRMS;
REAL8TimeSeries dataDouble;
REAL4TimeSeries dataSingle;

/* FUNCTION PROTOTYPES */
int WindowData(struct CommandLineArgsTag CLA);
int WindowDataTukey2(struct CommandLineArgsTag CLA);
int WindowDataHann(struct CommandLineArgsTag CLA);

int main(void) {

  const size_t NWINDOWS = 4;
  const size_t WINDOWLENGTH = 256 * 1800;

  // default command line arguments from MakeSFTs.c
  struct CommandLineArgsTag CLA = {
    .useSingle = 0,
    .windowR = 0.001,
  };

  // allocate memory
  char windownames[NWINDOWS][1024];
  REAL8Vector* windows[NWINDOWS];
  for (size_t i = 0; i < NWINDOWS; ++i) {
    windows[i] = XLALCreateREAL8Vector(WINDOWLENGTH);
    XLAL_CHECK_MAIN(windows[i] != NULL, XLAL_ENOMEM);
    for (size_t j = 0; j < WINDOWLENGTH; ++j) {
      windows[i]->data[j] = 1.0;
    }
  }

  size_t w = 0;

  {
    snprintf(windownames[w], sizeof(windownames[w]), "lalpulsar_MakeSFTs Matlab style Tukey window [windowR=%g]", CLA.windowR);
    dataDouble.data = windows[w];
    WindowData(CLA);
  }

  ++w;

  {
    snprintf(windownames[w], sizeof(windownames[w]), "lalpulsar_MakeSFTs Hann window");
    dataDouble.data = windows[w];
    WindowDataHann(CLA);
  }

  ++w;

  {
    REAL8 beta = CLA.windowR;
    REAL8Window *win = XLALCreateTukeyREAL8Window(WINDOWLENGTH, beta);
    XLAL_CHECK_MAIN(win != NULL, XLAL_EFUNC);
    snprintf(windownames[w], sizeof(windownames[w]), "XLALCreateTukeyREAL8Window(beta=%g)", beta);
    for (size_t j = 0; j < WINDOWLENGTH; ++j) {
      windows[w]->data[j] *= win->data->data[j];
    }
    XLALDestroyREAL8Window(win);
  }

  ++w;

  {
    REAL8Window *win = XLALCreateHannREAL8Window(WINDOWLENGTH);
    XLAL_CHECK_MAIN(win != NULL, XLAL_EFUNC);
    snprintf(windownames[w], sizeof(windownames[w]), "XLALCreateHannREAL8Window()");
    for (size_t j = 0; j < WINDOWLENGTH; ++j) {
      windows[w]->data[j] *= win->data->data[j];
    }
    XLALDestroyREAL8Window(win);
  }

  XLAL_CHECK_MAIN(w + 1 == NWINDOWS, XLAL_EFAILED);

  // output windows
  for (size_t i = 0; i < NWINDOWS; ++i) {
    printf("%s%c", windownames[i], i + 1 < NWINDOWS ? ',' : '\n');
  }
  for (size_t j = 0; j < WINDOWLENGTH; ++j) {
    for (size_t i = 0; i < NWINDOWS; ++i) {
      printf("%0.8f%c", windows[i]->data[j], i + 1 < NWINDOWS ? ',' : '\n');
    }
  }

  // cleanup
  for (size_t i = 0; i < NWINDOWS; ++i) {
    XLALDestroyREAL8Vector(windows[i]);
  }
  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

}
