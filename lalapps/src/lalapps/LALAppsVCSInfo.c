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
const char *const lalAppsVCSVersion = LALAPPS_VERSION;
const char *const lalAppsVCSId = LALAPPS_VCS_ID;
const char *const lalAppsVCSDate = LALAPPS_VCS_DATE;
const char *const lalAppsVCSBranch = LALAPPS_VCS_BRANCH;
const char *const lalAppsVCSTag = LALAPPS_VCS_TAG;
const char *const lalAppsVCSAuthor = LALAPPS_VCS_AUTHOR;
const char *const lalAppsVCSCommitter = LALAPPS_VCS_COMMITTER;
const char *const lalAppsVCSStatus = LALAPPS_VCS_STATUS;

/* global variables for vcs information - identable */
const char *const lalAppsVCSIdentId = LALAPPS_VCS_IDENT_ID;
const char *const lalAppsVCSIdentDate = LALAPPS_VCS_IDENT_DATE;
const char *const lalAppsVCSIdentBranch = LALAPPS_VCS_IDENT_BRANCH;
const char *const lalAppsVCSIdentTag = LALAPPS_VCS_IDENT_TAG;
const char *const lalAppsVCSIdentAuthor = LALAPPS_VCS_IDENT_AUTHOR;
const char *const lalAppsVCSIdentCommitter = LALAPPS_VCS_IDENT_COMMITTER;
const char *const lalAppsVCSIdentStatus = LALAPPS_VCS_IDENT_STATUS;

/* vcs information structure */
const struct tagLALVCSInfo lalAppsVCSInfo = { \
  LALAPPS_VCS_NAME, \
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
