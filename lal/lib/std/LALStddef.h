/*
*  Copyright (C) 2014, 2020 Karl Wette
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#ifndef _LALSTDDEF_H
#define _LALSTDDEF_H

/* forbid certain definitions in LAL code */
#if defined(LAL_STRICT_DEFS_ENABLED)
# if defined(__GNUC__)
/* assert() is often misused for error-checking code
   which should never be disabled in production usage */
#  pragma GCC poison assert
# endif
#endif

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

/* macros for compiler-specific attributes */
#if defined(__GNUC__)
# define _LAL_GCC_PRINTF_FORMAT_(NFMT, NARG) __attribute__ ((format (printf, NFMT, NARG)))
# define _LAL_GCC_VPRINTF_FORMAT_(NFMT) __attribute__ ((format (printf, NFMT, 0)))
#else
# define _LAL_GCC_PRINTF_FORMAT_(NFMT, NARG)
# define _LAL_GCC_VPRINTF_FORMAT_(NFMT)
#endif

#endif /* _LALSTDDEF_H */
