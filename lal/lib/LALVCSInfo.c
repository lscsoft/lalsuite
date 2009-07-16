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
 * Copyright (C) 2009 Adam Mercer
 */

#include <LALVCSInfo.h>

/* global variables for vcs information */
const char *lalVCSId = LAL_VCS_ID;
const char *lalVCSDate = LAL_VCS_DATE;
const char *lalVCSBranch = LAL_VCS_BRANCH;
const char *lalVCSTag = LAL_VCS_TAG;
const char *lalVCSAuthor = LAL_VCS_AUTHOR;
const char *lalVCSCommitter = LAL_VCS_COMMITTER;
const char *lalVCSStatus = LAL_VCS_STATUS;

/* library vcs information structure */
struct tagLALVCSInfo lalLibraryVCSInfo = {LAL_VCS_ID, \
  LAL_VCS_DATE, LAL_VCS_BRANCH, LAL_VCS_TAG, LAL_VCS_AUTHOR, \
  LAL_VCS_COMMITTER, LAL_VCS_STATUS};

/*
 * vim: tw=0 ts=2 et
 */
