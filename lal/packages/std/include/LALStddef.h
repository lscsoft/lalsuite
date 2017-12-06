/*
*  Copyright (C) 2013 Jolien Creighton
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

#ifndef _LALSTDDEF_H
#define _LALSTDDEF_H

/* macros for certain keywords */
#if __STDC_VERSION__ >= 199901L
# define _LAL_RESTRICT_ restrict
# define _LAL_INLINE_ inline
#elif defined __GNUC__
# define _LAL_RESTRICT_ __restrict__
# define _LAL_INLINE_ __inline__
#else
# define _LAL_RESTRICT_
# define _LAL_INLINE_
#endif

#endif /* _LALSTDDEF_H */
