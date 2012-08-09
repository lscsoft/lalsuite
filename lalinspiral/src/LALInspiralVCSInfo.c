/*
 * LALInspiralVCSInfo.c - LALInspiral VCS Information
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

#include <lal/LALInspiralConfig.h>
#include <lal/LALInspiralVCSInfo.h>

/* global variables for vcs information */
const char *const lalInspiralVCSVersion = LALINSPIRAL_VERSION;
const char *const lalInspiralVCSId = LALINSPIRAL_VCS_ID;
const char *const lalInspiralVCSDate = LALINSPIRAL_VCS_DATE;
const char *const lalInspiralVCSBranch = LALINSPIRAL_VCS_BRANCH;
const char *const lalInspiralVCSTag = LALINSPIRAL_VCS_TAG;
const char *const lalInspiralVCSAuthor = LALINSPIRAL_VCS_AUTHOR;
const char *const lalInspiralVCSCommitter = LALINSPIRAL_VCS_COMMITTER;
const char *const lalInspiralVCSStatus = LALINSPIRAL_VCS_STATUS;

/* global variables for vcs information - identable */
const char *const lalInspiralVCSIdentId = LALINSPIRAL_VCS_IDENT_ID;
const char *const lalInspiralVCSIdentDate = LALINSPIRAL_VCS_IDENT_DATE;
const char *const lalInspiralVCSIdentBranch = LALINSPIRAL_VCS_IDENT_BRANCH;
const char *const lalInspiralVCSIdentTag = LALINSPIRAL_VCS_IDENT_TAG;
const char *const lalInspiralVCSIdentAuthor = LALINSPIRAL_VCS_IDENT_AUTHOR;
const char *const lalInspiralVCSIdentCommitter = LALINSPIRAL_VCS_IDENT_COMMITTER;
const char *const lalInspiralVCSIdentStatus = LALINSPIRAL_VCS_IDENT_STATUS;

/* vcs information structure */
const struct tagLALVCSInfo lalInspiralVCSInfo = { \
  LALINSPIRAL_VCS_NAME, \
  LALINSPIRAL_VERSION, \
  LALINSPIRAL_VCS_ID, \
  LALINSPIRAL_VCS_DATE, \
  LALINSPIRAL_VCS_BRANCH, \
  LALINSPIRAL_VCS_TAG, \
  LALINSPIRAL_VCS_AUTHOR, \
  LALINSPIRAL_VCS_COMMITTER, \
  LALINSPIRAL_VCS_STATUS \
};

/*
 * vim: tw=0 ts=2 et
 */
