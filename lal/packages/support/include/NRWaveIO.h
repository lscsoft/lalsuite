/*
 * Copyright (C) 2006 S.Fairhurst, B. Krishnan, L.Santamaria
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

/** \defgroup NRWaveIO
 * \ingroup support
 * \author S.Fairhurst, B. Krishnan, L.Santamaria
 * 
 * \brief Module for reading/writing Numrel waveforms
 *

 *
 */
 
/** \file NRWaveIO.h
 *  \ingroup NRWaveIO   
 * \date $Date$
 *
 * 
 */

#ifndef _NRWAVEIO_H  	/* Double-include protection. */
#define _NRWAVEIO_H

/* includes */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>


#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID (NRWAVEIOH, "$Id$");


/** \name Error codes */
/*@{*/
#define NRWAVEIO_ENULL 	1
#define NRWAVEIO_EFILE 	2
#define NRWAVEIO_EHEADER 	3
#define NRWAVEIO_EVERSION 	4
#define NRWAVEIO_EVAL 		5
/*@}*/


REAL4TimeVectorSeries * XLALReadNRWave( REAL4  mass, CHAR   *filename); 


#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _SFTBIN_H */
