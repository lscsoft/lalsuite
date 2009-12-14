/*
*  Copyright (C) 2007 Cristina Valeria Torres, Jolien Creighton
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
 * File Name: TimeFreq.h
 *
 * New Maintainer: Torres, C (Univ TX at Browsville)
 * Author: E.C. Mottin
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * TimeFreq.h
 *
 * SYNOPSIS
 * #include <lal/TimeFreq.h>
 *
 * DESCRIPTION
 * Header file for the TFR package (computation of time-frequency
 * representation for the detection of gravitational waves from
 * unmodeled astrophysical sources)
 *
 * DIAGNOSTICS
 * ??
 *----------------------------------------------------------------------- */

 /*
 * 2. include-loop protection (see below). Note the naming convention!
 */

#ifndef _TIMEFREQ_H
#define _TIMEFREQ_H

/*
 * 3. Includes. This header may include others; if so, they go immediately
 *    after include-loop protection. Includes should appear in the following
 *    order:
 *    a. Standard library includes
 *    b. LDAS includes
 *    c. LAL includes
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (TIMEFREQH, "$Id$");

/*
 * 5. Macros. But, note that macros are deprecated.
 */

/*
 * 8. Structure, enum, union, etc., typdefs.
 */

#define CREATETFR_ENULL 1
#define CREATETFR_ENNUL 2
#define CREATETFR_EFROW 4
#define CREATETFR_ETCOL 8
#define CREATETFR_EMALL 16

#define CREATETFR_MSGENULL "Null pointer" /* ENULL */
#define CREATETFR_MSGENNUL "Non-null pointer" /* ENNUL */
#define CREATETFR_MSGEFROW "Illegal number of freq bins" /* EFROW */
#define CREATETFR_MSGETCOL "Illegal number of time instants" /* ETCOL */
#define CREATETFR_MSGEMALL "Malloc failure" /* EMALL */

#define DESTROYTFR_ENULL 1

#define DESTROYTFR_MSGENULL "Null pointer" /* ENULL */

#define CREATETFP_ENULL 1
#define CREATETFP_ENNUL 2
#define CREATETFP_EMALL 4
#define CREATETFP_EWSIZ 8
#define CREATETFP_ETYPE 16

#define CREATETFP_MSGENULL "Null pointer" /* ENULL */
#define CREATETFP_MSGENNUL "Non-null pointer" /* ENNUL */
#define CREATETFP_MSGEMALL "Malloc failure" /* EMALLOC */
#define CREATETFP_MSGEWSIZ "Invalid window length" /* ESIZ */
#define CREATETFP_MSGETYPE "Unknown TFR type" /* ETYPE */

#define DESTROYTFP_ENULL 1
#define DESTROYTFP_ETYPE 2

#define DESTROYTFP_MSGENULL "Null pointer" /* ENULL */
#define DESTROYTFP_MSGETYPE "Unknown TFR type" /* ETYPE */

#define TFR_ENULL 1
#define TFR_ENAME 2
#define TFR_EFROW 4
#define TFR_EWSIZ 8
#define TFR_ESAME 16
#define TFR_EBADT 32

#define TFR_MSGENULL "Null pointer" /* ENULL */
#define TFR_MSGENAME "TFR type mismatched" /* ENAME */
#define TFR_MSGEFROW "Invalid number of freq bins" /* EFROW */
#define TFR_MSGEWSIZ "Invalid window length" /* EWSIZ */
#define TFR_MSGESAME "Input/Output data vectors are the same" /* ESAME */
#define TFR_MSGEBADT "Invalid time instant" /* ETBAD */

/* Available TFR types */

#define TIME_FREQ_REP_NAMELIST {"Undefined","Spectrogram","WignerVille", "PSWignerVille","RSpectrogram"}

typedef enum tagTimeFreqRepType {
Undefined, Spectrogram, WignerVille, PSWignerVille, RSpectrogram
} TimeFreqRepType;

/* Time-Frequency Representation structure */

typedef struct tagTimeFreqRep {
  TimeFreqRepType type;             /* type of the TFR */
  INT4 fRow;                        /* number of freq bins in the TFR matrix */
  INT4 tCol;                        /* number of time bins in the TFR matrix */
  REAL4 *freqBin;	            /* freqs for each row of the matrix */
  INT4 *timeInstant;                /* time instants for each column of the TFR */
  REAL4 **map;                      /* TFR */
} TimeFreqRep;

/* TFR parameter structure */

typedef struct tagTimeFreqParam {
  TimeFreqRepType type;                   /* type of the TFR */
  REAL4Vector *windowT;                   /* (Sp, Rsp and Pswv) Window */
  REAL4Vector *windowF;                   /* (Pswv) Window */
} TimeFreqParam;

/* For memory allocation of the TFR and parameter structure */

typedef struct tagCreateTimeFreqIn {
  TimeFreqRepType type;             /* type of the TFR */
  INT4 fRow;                        /* number of freq bins in the TFR matrix */
  INT4 tCol;                        /* number of time bins in the TFR matrix */
  INT4 wlengthT;                    /* (Sp, Pswv and Rsp) Window length */
  INT4 wlengthF;                    /* (Pswv) Window */
} CreateTimeFreqIn;

/*
 * 9. Functions Declarations (i.e., prototypes).
 */

void LALCreateTimeFreqRep (LALStatus*, TimeFreqRep**, CreateTimeFreqIn*);
void LALCreateTimeFreqParam (LALStatus*, TimeFreqParam**, CreateTimeFreqIn*);
void LALDestroyTimeFreqRep (LALStatus*, TimeFreqRep**);
void LALDestroyTimeFreqParam (LALStatus*, TimeFreqParam**);
void LALTfrSp (LALStatus*, REAL4Vector*, TimeFreqRep*, TimeFreqParam*);
void LALTfrWv (LALStatus*, REAL4Vector*, TimeFreqRep*, TimeFreqParam*);
void LALTfrPswv (LALStatus*, REAL4Vector*, TimeFreqRep*, TimeFreqParam*);
void LALTfrRsp (LALStatus*, REAL4Vector*, TimeFreqRep*, TimeFreqParam*);
void LALDwindow (LALStatus*, REAL4Vector*, REAL4Vector*);

#ifdef  __cplusplus
}
#endif

#endif /* _TIMEFREQ_H */
