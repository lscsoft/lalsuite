/*
 * LALBurstVCSInfo.c - LALBurst VCS Information
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

#include <lal/LALBurstConfig.h>
#include <lal/LALBurstVCSInfo.h>

/* global variables for vcs information */
const char *lalBurstVCSVersion = LALBURST_VERSION;
const char *lalBurstVCSId = LALBURST_VCS_ID;
const char *lalBurstVCSDate = LALBURST_VCS_DATE;
const char *lalBurstVCSBranch = LALBURST_VCS_BRANCH;
const char *lalBurstVCSTag = LALBURST_VCS_TAG;
const char *lalBurstVCSAuthor = LALBURST_VCS_AUTHOR;
const char *lalBurstVCSCommitter = LALBURST_VCS_COMMITTER;
const char *lalBurstVCSStatus = LALBURST_VCS_STATUS;

/* global variables for vcs information - identable */
const char *lalBurstVCSIdentId = LALBURST_VCS_IDENT_ID;
const char *lalBurstVCSIdentDate = LALBURST_VCS_IDENT_DATE;
const char *lalBurstVCSIdentBranch = LALBURST_VCS_IDENT_BRANCH;
const char *lalBurstVCSIdentTag = LALBURST_VCS_IDENT_TAG;
const char *lalBurstVCSIdentAuthor = LALBURST_VCS_IDENT_AUTHOR;
const char *lalBurstVCSIdentCommitter = LALBURST_VCS_IDENT_COMMITTER;
const char *lalBurstVCSIdentStatus = LALBURST_VCS_IDENT_STATUS;

/* vcs information structure */
struct tagLALVCSInfo lalBurstVCSInfo = { \
  LALBURST_VCS_NAME, \
  LALBURST_VERSION, \
  LALBURST_VCS_ID, \
  LALBURST_VCS_DATE, \
  LALBURST_VCS_BRANCH, \
  LALBURST_VCS_TAG, \
  LALBURST_VCS_AUTHOR, \
  LALBURST_VCS_COMMITTER, \
  LALBURST_VCS_STATUS \
};

/*
 * vim: tw=0 ts=2 et
 */
