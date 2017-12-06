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
#ifndef _FOLDAMPLITUDES_H
#define _FOLDAMPLITUDES_H

#include <lal/LALStdlib.h>
/* ****** INCLUDE ANY OTHER LAL HEADERS needed for header (NOT module) *** */

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \author Mendell, Greg A.
 * \defgroup FoldAmplitudes_h Header FoldAmplitudes.h
 * \ingroup pkg_pulsarFold
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/FoldAmplitudes.h>
 * \endcode
 *
 * Contains prototypes for:
 *
 * struct FoldAmplitudesInput
 *
 * struct FoldAmplitudesParams
 *
 * function LALFoldAmplitudes
 *
 * %[Generic documentation on the header; this is the main place to
 * %document any stuff not specific to the module]
 *
 * ### Structures ###
 *
 * %[Document here any structures defined in the header.
 * %Also include any of them in the index; e.g.:]
 * \code
 *
 * typedef struct tagFoldAmplitudesInput
 * {
 *
 * REAL4Vector  *amplitudeVec;  input vector of amplitudes
 * REAL4Vector          *phaseVec;      input vector of phases
 *
 * } FoldAmplitudesInput;
 *
 * typedef struct tagFoldAmplitudesParams
 * {
 *
 * INT4		numBins;      number of bins
 * REAL4		binMin;       minimum phase to bin
 * REAL4		binMax;       maximum phase to bin
 *
 * } FoldAmplitudesParams;
 *
 * \endcode
 *
 */
/*@{*/

/**\name Error Codes */
/*@{*/
#define FOLDAMPLITUDESH_ENULLP          1
#define FOLDAMPLITUDESH_EVECSIZE        2
#define FOLDAMPLITUDESH_ENUMBINS        3
#define FOLDAMPLITUDESH_EBINSIZE        4
#define FOLDAMPLITUDESH_EBINMIN         5

#define FOLDAMPLITUDESH_MSGENULLP       "Null pointer!"
#define FOLDAMPLITUDESH_MSGEVECSIZE     "Input vectors were not the same length!"
#define FOLDAMPLITUDESH_MSGENUMBINS     "Number of bins was less than 1!"
#define FOLDAMPLITUDESH_MSGEBINSIZE     "Bin max was less than bin min!"
#define FOLDAMPLITUDESH_MSGEBINMIN      "Bin min was not zero; nonzero bin min has not yet been implemented!"
/*@}*/

/* ***** DEFINE OTHER GLOBAL CONSTANTS OR MACROS *********** */

/* ***** DEFINE NEW STRUCTURES AND TYPES *********** */

typedef struct tagFoldAmplitudesInput
{

  REAL4Vector   *amplitudeVec;  /* input vector of amplitudes */
  REAL4Vector   *phaseVec;      /* input vector of phases */

} FoldAmplitudesInput;

typedef struct tagFoldAmplitudesParams
{

  INT4		numBins;       /* number of bins */
  REAL4		binMin;        /* minimum phase to bin */
  REAL4		binMax;        /* maximum phase to bin */

} FoldAmplitudesParams;

/****** INCLUDE EXTERNAL GLOBAL VARIABLES ************/

/*
LALFoldAmplitudes:
inputs: a vector of amplitudes and a vector of phases.
params: number of phase bins,  the minimum phase to bin, and the maximum phase to bin.
action: for each phase, the phase is first reduced by modulo arithmetic to a value
        between binMin and binMax.  The corresponding amplitudes is then added to
        the corresponding phase bins.  The width of each bin is (binMax - binMin)/numBins.
output: a vector of folded amplitudes; component i is the folded amplitude for phase bin i.
*/

void
LALFoldAmplitudes( LALStatus                      *status,
                     REAL4Vector                  *output,
                     const FoldAmplitudesInput    *input,
                     const FoldAmplitudesParams   *params );

/*@}*/

#ifdef  __cplusplus
}
#endif

#endif /* _FOLDAMPLITUDES_H */
