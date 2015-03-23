/*
 * Copyright (C) 2005,2007,2008,2010,2015  Kipp Cannon
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


#ifndef _THRESHOLDS_H
#define _THRESHOLDS_H


#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

/**
 * \defgroup Thresholds_h Header Thresholds.h
 * \ingroup lalburst_burstsearch
 *
 * \brief UNDOCUMENTED
 */
  /*@{*/


double XLALLogChisqCCDF(
	double chi2,
	double dof
);


  /*@}*/

#ifdef  __cplusplus
}
#endif  /* C++ protection. */


#endif
