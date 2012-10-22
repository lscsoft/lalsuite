/*
 * LALDetCharVCSInfo.c - LALDetChar VCS Information
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

#include <lal/LALDetCharConfig.h>
#include <lal/LALDetCharVCSInfo.h>

/* global variables for vcs information */
const char *const lalDetCharVCSVersion = LALDETCHAR_VERSION;
const char *const lalDetCharVCSId = LALDETCHAR_VCS_ID;
const char *const lalDetCharVCSDate = LALDETCHAR_VCS_DATE;
const char *const lalDetCharVCSBranch = LALDETCHAR_VCS_BRANCH;
const char *const lalDetCharVCSTag = LALDETCHAR_VCS_TAG;
const char *const lalDetCharVCSAuthor = LALDETCHAR_VCS_AUTHOR;
const char *const lalDetCharVCSCommitter = LALDETCHAR_VCS_COMMITTER;
const char *const lalDetCharVCSStatus = LALDETCHAR_VCS_STATUS;

/* global variables for vcs information - identable */
const char *const lalDetCharVCSIdentId = LALDETCHAR_VCS_IDENT_ID;
const char *const lalDetCharVCSIdentDate = LALDETCHAR_VCS_IDENT_DATE;
const char *const lalDetCharVCSIdentBranch = LALDETCHAR_VCS_IDENT_BRANCH;
const char *const lalDetCharVCSIdentTag = LALDETCHAR_VCS_IDENT_TAG;
const char *const lalDetCharVCSIdentAuthor = LALDETCHAR_VCS_IDENT_AUTHOR;
const char *const lalDetCharVCSIdentCommitter = LALDETCHAR_VCS_IDENT_COMMITTER;
const char *const lalDetCharVCSIdentStatus = LALDETCHAR_VCS_IDENT_STATUS;

/* vcs information structure */
const struct tagLALVCSInfo lalDetCharVCSInfo = { \
  LALDETCHAR_VCS_NAME, \
  LALDETCHAR_VERSION, \
  LALDETCHAR_VCS_ID, \
  LALDETCHAR_VCS_DATE, \
  LALDETCHAR_VCS_BRANCH, \
  LALDETCHAR_VCS_TAG, \
  LALDETCHAR_VCS_AUTHOR, \
  LALDETCHAR_VCS_COMMITTER, \
  LALDETCHAR_VCS_STATUS \
};

/*
 * vim: tw=0 ts=2 et
 */
