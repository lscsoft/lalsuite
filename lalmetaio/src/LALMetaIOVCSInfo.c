/*
 * LALMetaIOVCSInfo.c - LALMetaIO VCS Information
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

#include <lal/LALVCSInfo.h>

#include <lal/LALMetaIOConfig.h>
#include <lal/LALMetaIOVCSInfo.h>

/* global variables for vcs information */
const char *const lalMetaIOVCSVersion = LALMETAIO_VERSION;
const char *const lalMetaIOVCSId = LALMETAIO_VCS_ID;
const char *const lalMetaIOVCSDate = LALMETAIO_VCS_DATE;
const char *const lalMetaIOVCSBranch = LALMETAIO_VCS_BRANCH;
const char *const lalMetaIOVCSTag = LALMETAIO_VCS_TAG;
const char *const lalMetaIOVCSAuthor = LALMETAIO_VCS_AUTHOR;
const char *const lalMetaIOVCSCommitter = LALMETAIO_VCS_COMMITTER;
const char *const lalMetaIOVCSStatus = LALMETAIO_VCS_STATUS;

/* global variables for vcs information - identable */
const char *const lalMetaIOVCSIdentId = LALMETAIO_VCS_IDENT_ID;
const char *const lalMetaIOVCSIdentDate = LALMETAIO_VCS_IDENT_DATE;
const char *const lalMetaIOVCSIdentBranch = LALMETAIO_VCS_IDENT_BRANCH;
const char *const lalMetaIOVCSIdentTag = LALMETAIO_VCS_IDENT_TAG;
const char *const lalMetaIOVCSIdentAuthor = LALMETAIO_VCS_IDENT_AUTHOR;
const char *const lalMetaIOVCSIdentCommitter = LALMETAIO_VCS_IDENT_COMMITTER;
const char *const lalMetaIOVCSIdentStatus = LALMETAIO_VCS_IDENT_STATUS;

/* vcs information structure */
const struct tagLALVCSInfo lalMetaIOVCSInfo = { \
  LALMETAIO_VCS_NAME, \
  LALMETAIO_VERSION, \
  LALMETAIO_VCS_ID, \
  LALMETAIO_VCS_DATE, \
  LALMETAIO_VCS_BRANCH, \
  LALMETAIO_VCS_TAG, \
  LALMETAIO_VCS_AUTHOR, \
  LALMETAIO_VCS_COMMITTER, \
  LALMETAIO_VCS_STATUS \
};

/*
 * vim: tw=0 ts=2 et
 */
