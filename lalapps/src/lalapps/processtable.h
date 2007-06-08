/*
*  Copyright (C) 2007 Duncan Brown, Kipp Cannon
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

#ifndef PROCESSTABLE_H_
#define PROCESSTABLE_H_

#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>

#ifdef  __cplusplus
extern "C" {
#endif

void
populate_process_table (
    LALStatus           *status,
    ProcessTable        *ptable,
    const CHAR          *program,
    const CHAR          *cvs_revision,
    const CHAR          *cvs_source,
    const CHAR          *cvs_date
    );

#ifdef  __cplusplus
}
#endif

#endif
