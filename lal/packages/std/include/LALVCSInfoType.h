/*
 * LALVCSInfoType.h - LAL VCS Information Type
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
 */

#ifndef _LALVCSINFOTYPE_H
#define _LALVCSINFOTYPE_H

#ifdef __cplusplus
extern "C" {
#endif

/* define vcs information structure */
typedef struct tagLALVCSInfo
{
  const char *const name;
  const char *const version;
  const char *const vcsId;
  const char *const vcsDate;
  const char *const vcsBranch;
  const char *const vcsTag;
  const char *const vcsAuthor;
  const char *const vcsCommitter;
  const char *const vcsStatus;
} LALVCSInfo;

/* function to compare two LALVCSInfo structures */
int XLALVCSInfoCompare(const LALVCSInfo *header, const LALVCSInfo *library);

#ifdef __cplusplus
}
#endif

#endif /* _LALVCSINFOTYPE_H */

/*
 * vim: tw=0 ts=2 et
 */
