/*
 * LALXMLLibVCSInfo.h - LALXML VCS Information Header
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

#ifndef _LALXMLLIBVCSINFO_H
#define _LALXMLLIBVCSINFO_H

#include <lal/LALLibVCSInfo.h>
#include <lal/LALXMLConfig.h>

#ifdef __cplusplus
extern "C" {
#endif

/* global variables for vcs information, defined in LALXMLVCSInfo.c */
extern const char *const lalXMLVCSVersion;
extern const char *const lalXMLVCSId;
extern const char *const lalXMLVCSDate;
extern const char *const lalXMLVCSBranch;
extern const char *const lalXMLVCSTag;
extern const char *const lalXMLVCSAuthor;
extern const char *const lalXMLVCSCommitter;
extern const char *const lalXMLVCSStatus;

/* vcs information structures */
extern const struct tagLALVCSInfo lalXMLVCSInfo;

#ifdef __cplusplus
}
#endif

#endif /* _LALXMLLIBVCSINFO_H */

/*
 * vim: tw=0 ts=2 et
 */
