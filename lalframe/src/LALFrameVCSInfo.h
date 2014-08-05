/*
 * LALFrameVCSInfo.h - LALFrame VCS Information
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
 * Copyright (C) 2009-2013 Adam Mercer
 * Copyright (C) 2014 Karl Wette
 */

#ifndef _LALFRAMEVCSINFO_H
#define _LALFRAMEVCSINFO_H

#include <lal/LALVCSInfoType.h>

#ifdef __cplusplus
extern "C" {
#endif

/* global variables for vcs information */
extern const char *const lalFrameVCSVersion;
extern const char *const lalFrameVCSId;
extern const char *const lalFrameVCSDate;
extern const char *const lalFrameVCSBranch;
extern const char *const lalFrameVCSTag;
extern const char *const lalFrameVCSAuthor;
extern const char *const lalFrameVCSCommitter;
extern const char *const lalFrameVCSStatus;

/* global variables for vcs information - identable */
extern const char *const lalFrameVCSIdentId;
extern const char *const lalFrameVCSIdentDate;
extern const char *const lalFrameVCSIdentBranch;
extern const char *const lalFrameVCSIdentTag;
extern const char *const lalFrameVCSIdentAuthor;
extern const char *const lalFrameVCSIdentCommitter;
extern const char *const lalFrameVCSIdentStatus;

/* library vcs information structure */
extern const struct tagLALVCSInfo lalFrameVCSInfo;

/* configure arguments */
extern const char *const lalFrameConfigureArgs;

/* configure date */
extern const char *const lalFrameConfigureDate;

/* build date */
extern const char *const lalFrameBuildDate;

#ifdef __cplusplus
}
#endif

#endif /* _LALFRAMEVCSINFO_H */

/*
 * vim: tw=0 ts=2 et
 */
