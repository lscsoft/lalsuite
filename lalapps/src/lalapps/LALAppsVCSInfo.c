/*
 * LALAppsVCSInfo.c - LALApps VCS Information
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

#include <config.h>

#include <lal/LALVCSInfo.h>

#include <LALAppsVCSInfo.h>

/* global variables for vcs information */
const char *lalAppsVCSVersion = LALAPPS_VERSION;
const char *lalAppsVCSId = LALAPPS_VCS_ID;
const char *lalAppsVCSDate = LALAPPS_VCS_DATE;
const char *lalAppsVCSBranch = LALAPPS_VCS_BRANCH;
const char *lalAppsVCSTag = LALAPPS_VCS_TAG;
const char *lalAppsVCSAuthor = LALAPPS_VCS_AUTHOR;
const char *lalAppsVCSCommitter = LALAPPS_VCS_COMMITTER;
const char *lalAppsVCSStatus = LALAPPS_VCS_STATUS;

/* global variables for vcs information - identable */
const char *lalAppsVCSIdentId = LALAPPS_VCS_IDENT_ID;
const char *lalAppsVCSIdentDate = LALAPPS_VCS_IDENT_DATE;
const char *lalAppsVCSIdentBranch = LALAPPS_VCS_IDENT_BRANCH;
const char *lalAppsVCSIdentTag = LALAPPS_VCS_IDENT_TAG;
const char *lalAppsVCSIdentAuthor = LALAPPS_VCS_IDENT_AUTHOR;
const char *lalAppsVCSIdentCommitter = LALAPPS_VCS_IDENT_COMMITTER;
const char *lalAppsVCSIdentStatus = LALAPPS_VCS_IDENT_STATUS;

/* vcs information structure */
struct tagLALVCSInfo lalAppsVCSInfo = { \
  LALAPPS_VERSION, \
  LALAPPS_VCS_ID, \
  LALAPPS_VCS_DATE, \
  LALAPPS_VCS_BRANCH, \
  LALAPPS_VCS_TAG, \
  LALAPPS_VCS_AUTHOR, \
  LALAPPS_VCS_COMMITTER, \
  LALAPPS_VCS_STATUS \
};

/*
 * vim: tw=0 ts=2 et
 */
