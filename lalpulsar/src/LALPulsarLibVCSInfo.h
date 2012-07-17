/*
 * LALPulsarLibVCSInfo.h - LALPulsar VCS Information Header
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

#ifndef _LALPULSARLIBVCSINFO_H
#define _LALPULSARLIBVCSINFO_H

#include <lal/LALLibVCSInfo.h>
#include <lal/LALPulsarConfig.h>

#ifdef __cplusplus
extern "C" {
#endif

/* global variables for vcs information, defined in LALPulsarVCSInfo.c */
extern const char *const lalPulsarVCSVersion;
extern const char *const lalPulsarVCSId;
extern const char *const lalPulsarVCSDate;
extern const char *const lalPulsarVCSBranch;
extern const char *const lalPulsarVCSTag;
extern const char *const lalPulsarVCSAuthor;
extern const char *const lalPulsarVCSCommitter;
extern const char *const lalPulsarVCSStatus;

/* vcs information structures */
extern const struct tagLALVCSInfo lalPulsarVCSInfo;

#ifdef __cplusplus
}
#endif

#endif /* _LALPULSARLIBVCSINFO_H */

/*
 * vim: tw=0 ts=2 et
 */
