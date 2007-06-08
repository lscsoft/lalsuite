/*
*  Copyright (C) 2007 B.S. Sathyaprakash
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

/* <lalVerbatim file="FrameReadHV"> 

Author: Sathyaprakash, B.S.

</lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{GEO600FrameChannels.h}}
\label{s:GEO600FrameChannels.h}

Header file for frmae lib interface functions.

\subsection*{Synopsis}
\begin{verbatim}
#include "GEO600FrameChannels.h"
\end{verbatim}

\noindent This header file covers routines that are used in 
reading frame files.

</lalLaTeX> */

/*
  All functions return 0 if successful but negative numbers if unsuccessful
  They all print error messages to standard output;
*/


/* Sets base directory */
int
SetBaseDirectory(char *);
/*
#1 : the name of the base directory where the GEO Data resides
*/


/* gets Calibration Function */
int
GetCalibFunction(int,char *, float *, float, int);
/*
#1 : the gpsTime for which u want information : must be a multiple of 60
#2 : the name of the channel
#3 : the array in which u  want the information: user allocates memory (size)
#4 : the fundamental frequency
#5 : the size of the array
*/


/* Gets single Channel Data */

int
GetChannelData(int,int,char *, float *, int);
/* 
#1 : the gpsTime for which u want information  
#2 : the number of seconds of Data the user wants
#3 : the name of the channel
#4 : the array in which user  want the information:user allocates memory (size)
#5 : the size of the array allocated 
*/

/* Get double precision data  */
int
GetChannelDoubleData(int,int,char *, double *, int);
/* 
#1 : the gpsTime for which u want information  
#2 : the number of seconds of Data the user wants
#3 : the name of the channel
#4 : the array in which user  want the information:user allocates memory (size)
#5 : the size of the array allocated 
*/

/* gets Sample Rate */

int
GetSampleRate(int, char *, float *);
/*
#1 : the gpsTime for which u want information  
#2 : the name of the channel
#3 : a pointer to a single float containing the sampleRate 
*/

/* Gets White Noise */
int
GetWhiteNoiseData(int, float *);
/*
#1 : the number of points u want
#2 : a pointer to a array (user allocates memory 
*/ 
