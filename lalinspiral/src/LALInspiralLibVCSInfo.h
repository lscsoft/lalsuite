/*
 * LALInspiralLibVCSInfo.h - LALInspiral VCS Information Header
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
 * Copyright (C) 2009-2011 Adam Mercer
 */

#ifndef _LALINSPIRALLIBVCSINFO_H
#define _LALINSPIRALLIBVCSINFO_H

#include <lal/LALLibVCSInfo.h>
#include <lal/LALInspiralConfig.h>

#ifdef __cplusplus
extern "C" {
#endif

/* global variables for vcs information, defined in LALInspiralVCSInfo.c */
extern const char *const lalInspiralVCSVersion;
extern const char *const lalInspiralVCSId;
extern const char *const lalInspiralVCSDate;
extern const char *const lalInspiralVCSBranch;
extern const char *const lalInspiralVCSTag;
extern const char *const lalInspiralVCSAuthor;
extern const char *const lalInspiralVCSCommitter;
extern const char *const lalInspiralVCSStatus;

/* vcs information structures */
extern const struct tagLALVCSInfo lalInspiralVCSInfo;

#ifdef __cplusplus
}
#endif

#endif /* _LALINSPIRALLIBVCSINFO_H */

/*
 * vim: tw=0 ts=2 et
 */
