/*
 * LALVCSInfo.c - LAL VCS Information
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
 * Copyright (C) 2009,2010 Adam Mercer
 */

#include <string.h>
#include "config.h"
#include <LALVCSInfo.h>

/* global variables for vcs information */
const char *const lalVCSVersion = LAL_VERSION;
const char *const lalVCSId = LAL_VCS_ID;
const char *const lalVCSDate = LAL_VCS_DATE;
const char *const lalVCSBranch = LAL_VCS_BRANCH;
const char *const lalVCSTag = LAL_VCS_TAG;
const char *const lalVCSAuthor = LAL_VCS_AUTHOR;
const char *const lalVCSCommitter = LAL_VCS_COMMITTER;
const char *const lalVCSStatus = LAL_VCS_STATUS;

/* global variables for vcs information - identable */
const char *const lalVCSIdentId = LAL_VCS_IDENT_ID;
const char *const lalVCSIdentDate = LAL_VCS_IDENT_DATE;
const char *const lalVCSIdentBranch = LAL_VCS_IDENT_BRANCH;
const char *const lalVCSIdentTag = LAL_VCS_IDENT_TAG;
const char *const lalVCSIdentAuthor = LAL_VCS_IDENT_AUTHOR;
const char *const lalVCSIdentCommitter = LAL_VCS_IDENT_COMMITTER;
const char *const lalVCSIdentStatus = LAL_VCS_IDENT_STATUS;

/* library vcs information structure */
const struct tagLALVCSInfo lalVCSInfo = { \
  LAL_NAME, \
  LAL_VERSION, \
  LAL_VCS_ID, \
  LAL_VCS_DATE, \
  LAL_VCS_BRANCH, \
  LAL_VCS_TAG, \
  LAL_VCS_AUTHOR, \
  LAL_VCS_COMMITTER, \
  LAL_VCS_STATUS \
};

/* function to compare two LALVCSInfo structures */
int XLALVCSInfoCompare(const LALVCSInfo *header, const LALVCSInfo *library)
{
  /* check for header/library mismatch */
  if (strcmp(header->name, library->name) || \
      strcmp(header->version, library->version) || \
      strcmp(header->vcsId, library->vcsId) || \
      strcmp(header->vcsDate, library->vcsDate) || \
      strcmp(header->vcsBranch, library->vcsBranch) || \
      strcmp(header->vcsTag, library->vcsTag) || \
      strcmp(header->vcsAuthor, library->vcsAuthor) || \
      strcmp(header->vcsCommitter, library->vcsCommitter) || \
      strcmp(header->vcsStatus, library->vcsStatus))
  {
    /* version mismatch */
    return 1;
  }

  return 0;
}

/*
 * vim: tw=0 ts=2 et
 */
