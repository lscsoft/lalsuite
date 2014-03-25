/*
 * LALInspiralBuildInfo.h - LALInspiral Build Information Header
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA 02111-1307 USA
 *
 * Copyright (C) 2013 Karl Wette
 */

#ifndef _LALINSPIRALBUILDINFO_H
#define _LALINSPIRALBUILDINFO_H

#ifdef __cplusplus
extern "C" {
#endif

/* configure arguments */
extern const char *const lalInspiralConfigureArgs;

/* configure date */
extern const char *const lalInspiralConfigureDate;

/* build date */
extern const char *const lalInspiralBuildDate;

#ifdef __cplusplus
}
#endif

#endif /* _LALINSPIRALBUILDINFO_H */
