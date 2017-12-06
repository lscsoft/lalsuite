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

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Resample.h>

/**
 * \author Creighton, T. D.
 * \ingroup Resample_h
 * \brief Destroys an object of type ResampleRules, and sets <tt>*rules</tt> to \c NULL.
 */
void
LALDestroyResampleRules( LALStatus     *stat,
			 ResampleRules **rules )
{
  INITSTATUS(stat);

  /* Make sure that the handle is non-null, that it points to a
     non-null pointer, and that the interval and shift fields are
     non-null.) */
  ASSERT(rules,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(*rules,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT((*rules)->interval,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT((*rules)->shift,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);

  /* Free all the memory. */
  LALFree((*rules)->interval);
  LALFree((*rules)->shift);
  LALFree(*rules);

  /* Point the handle to NULL, then get out of here. */
  *rules=NULL;
  RETURN(stat);
}
