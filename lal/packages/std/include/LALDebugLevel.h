/*
*  Copyright (C) 2013 Jolien Creighton, Kipp Cannon
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

#ifndef _LALDEBUGLEVEL_H
#define _LALDEBUGLEVEL_H

#ifdef __cplusplus
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

/** lalDebugLevel bit field values */
enum {
    LALERRORBIT = 0001,   /**< enable error messages */
    LALWARNINGBIT = 0002, /**< enable warning messages */
    LALINFOBIT = 0004,    /**< enable info messages */
    LALTRACEBIT = 0010,   /**< enable tracing messages */
    LALMEMDBGBIT = 0020,  /**< enable memory debugging routines */
    LALMEMPADBIT = 0040,  /**< enable memory padding */
    LALMEMTRKBIT = 0100,  /**< enable memory tracking */
    LALMEMINFOBIT = 0200  /**< enable memory info messages */
};

/** composite lalDebugLevel values */
enum {
    LALNDEBUG = 0,      /**< no debug */
    LALERROR = LALERRORBIT,             /**< enable error messages */
    LALWARNING = LALWARNINGBIT,         /**< enable warning messages */
    LALINFO = LALINFOBIT,               /**< enable info messages */
    LALTRACE = LALTRACEBIT,             /**< enable tracing messages */
    LALMSGLVL1 = LALERRORBIT,           /**< enable error messages */
    LALMSGLVL2 = LALERRORBIT | LALWARNINGBIT,   /**< enable error and warning messages */
    LALMSGLVL3 = LALERRORBIT | LALWARNINGBIT | LALINFOBIT,      /**< enable error, warning, and info messages */
    LALMEMDBG = LALMEMDBGBIT | LALMEMPADBIT | LALMEMTRKBIT,     /**< enable memory debugging tools */
    LALMEMTRACE = LALTRACEBIT | LALMEMDBG | LALMEMINFOBIT,      /**< enable memory tracing tools */
    LALALLDBG = ~LALNDEBUG      /**< enable all debugging */
};

#define lalDebugLevel (XLALGetDebugLevel())
int XLALGetDebugLevel(void);

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif
#endif /* _LALDEBUGLEVEL_H */
