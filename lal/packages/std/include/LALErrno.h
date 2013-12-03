/*
*  Copyright (C) 2007 Bernd Machenschalk, Kipp Cannon
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
 * \author: Cannon, K. C.
 */

#ifndef _LALERRNO_H
#define _LALERRNO_H

/*
 * Error codes and corresponding error messages.
 */

#define LAL_FAIL_ERR	1
#define LAL_FAIL_MSG	"operation failed"
#define LAL_NULL_ERR	2
#define LAL_NULL_MSG	"unexpected NULL pointer"
#define LAL_NNULL_ERR	3
#define LAL_NNULL_MSG	"unexpected non-NULL pointer"
#define LAL_NOMEM_ERR	4
#define LAL_NOMEM_MSG	"out of memory"
#define LAL_RANGE_ERR	5
#define LAL_RANGE_MSG	"parameter out of range"
#define LAL_BADPARM_ERR 6
#define LAL_BADPARM_MSG "invalid parameter value"

#endif /* _LALERRNO_H */
