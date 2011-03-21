/*
 * $Id$
 *
 * Copyright (C) 2007  Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#ifndef _SEQUENCE_H
#define _SEQUENCE_H


#include <stddef.h>
#include <lal/LALDatatypes.h>
#include <lal/LALRCSID.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

NRCSID(SEQUENCEH, "$Id:");

define(`DATATYPE',COMPLEX8)
define(`SQUAREDATATYPE',REAL4)
include(SequenceH.m4)

define(`DATATYPE',COMPLEX16)
define(`SQUAREDATATYPE',REAL8)
include(SequenceH.m4)

define(`DATATYPE',REAL4)
define(`SQUAREDATATYPE',REAL4)
include(SequenceH.m4)

define(`DATATYPE',REAL8)
define(`SQUAREDATATYPE',REAL8)
include(SequenceH.m4)

define(`DATATYPE',INT2)
define(`SQUAREDATATYPE',UINT2)
include(SequenceH.m4)

define(`DATATYPE',INT4)
define(`SQUAREDATATYPE',UINT4)
include(SequenceH.m4)

define(`DATATYPE',INT8)
define(`SQUAREDATATYPE',UINT8)
include(SequenceH.m4)

define(`DATATYPE',UINT2)
define(`SQUAREDATATYPE',UINT2)
include(SequenceH.m4)

define(`DATATYPE',UINT4)
define(`SQUAREDATATYPE',UINT4)
include(SequenceH.m4)

define(`DATATYPE',UINT8)
define(`SQUAREDATATYPE',UINT8)
include(SequenceH.m4)

#ifdef __cplusplus
#pragma {
}
#endif

#endif  /* _SEQUENCE_H */
