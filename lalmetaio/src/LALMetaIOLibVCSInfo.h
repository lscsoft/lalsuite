/*
 * LALMetaIOLibVCSInfo.h - LALMetaIO VCS Information Header
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
 * Copyright (C) 2009-2010 Adam Mercer
 */

#ifndef _LALMETAIOLIBVCSINFO_H
#define _LALMETAIOLIBVCSINFO_H

#include <lal/LALLibVCSInfo.h>
#include <lal/LALMetaIOConfig.h>

#ifdef __cplusplus
extern "C" {
#endif

/* global variables for vcs information, defined in LALMetaIOVCSInfo.c */
extern const char *const lalMetaIOVCSVersion;
extern const char *const lalMetaIOVCSId;
extern const char *const lalMetaIOVCSDate;
extern const char *const lalMetaIOVCSBranch;
extern const char *const lalMetaIOVCSTag;
extern const char *const lalMetaIOVCSAuthor;
extern const char *const lalMetaIOVCSCommitter;
extern const char *const lalMetaIOVCSStatus;

/* vcs information structures */
extern const struct tagLALVCSInfo lalMetaIOVCSInfo;

#ifdef __cplusplus
}
#endif

#endif /* _LALMETAIOLIBVCSINFO_H */

/*
 * vim: tw=0 ts=2 et
 */
