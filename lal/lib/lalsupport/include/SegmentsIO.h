/*
*  Copyright (C) 2007 Jolien Creighton, Peter Shawhan
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

/**
 * \addtogroup SegmentsIO_h
 * \author Peter Shawhan
 *
 * \brief Provides segment list reading and writing functions as part of the \c lalsupport library.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/Segments.h>
 * #include <lal/SegmentsIO.h>
 * \endcode
 *
 * ### Notes ###
 *
 * The baseline format of a segment list file is described at
 * <tt>http://www.lsc-group.phys.uwm.edu/daswg/docs/technical/seglist_format.html</tt> .
 */
/*@{*/
/**\name Error Codes */ /*@{*/
#define SEGMENTSIOH_ENULL 1     /**< Null pointer passed to function */
#define SEGMENTSIOH_EINVAL 2    /**< LALSegList structure was not properly initialized */
#define SEGMENTSIOH_EBADOPT 3   /**< Invalid option letter in options string */
#define SEGMENTSIOH_ENOFMT 4    /**< No output format specified in options string */
#define SEGMENTSIOH_EOPENR 5    /**< Error opening segment list file for reading */
#define SEGMENTSIOH_EOPENW 6    /**< Error opening segment list file for writing */
#define SEGMENTSIOH_EFMT 7      /**< Segment list file is not in a recognized format */
#define SEGMENTSIOH_EPARSE 8    /**< Parsing error while reading from file */
#define SEGMENTSIOH_EDOM 9      /**< GPS times do not represent a valid segment */
#define SEGMENTSIOH_EINT 10     /**< Internal error in SegmentsIO module */
/*@}*/
/*@}*/

#define SEGMENTSIOH_MSGENULL "Null pointer passed to function"
#define SEGMENTSIOH_MSGEINVAL "LALSegList structure was not properly initialized"
#define SEGMENTSIOH_MSGEBADOPT "Invalid option letter in options string"
#define SEGMENTSIOH_MSGENOFMT "No output format specified in options string"
#define SEGMENTSIOH_MSGEOPENR "Error opening segment list file for reading"
#define SEGMENTSIOH_MSGEOPENW "Error opening segment list file for writing"
#define SEGMENTSIOH_MSGEFMT "Segment list file is not in a recognized format"
#define SEGMENTSIOH_MSGEPARSE "Parsing error while reading from file"
#define SEGMENTSIOH_MSGEDOM "GPS times do not represent a valid segment"
#define SEGMENTSIOH_MSGEINT "Internal error in SegmentsIO module"

#ifndef _SEGMENTSIO_H
#define _SEGMENTSIO_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <lal/FileIO.h>
#include <lal/Segments.h>

/* Function prototypes */

void
LALSegListRead( LALStatus *status, LALSegList *seglist, const CHAR *fileName, const CHAR *options );

void
LALSegListWrite( LALStatus *status, LALSegList *seglist, const CHAR *fileName, const CHAR *options );

#ifdef __cplusplus
}
#endif

#endif /* _SEGMENTSIO_H */
