/*
 * LALInferenceVCSInfo.c - LALInference VCS Information
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

#include <lal/LALInferenceConfig.h>
#include <lal/LALInferenceVCSInfo.h>

/* global variables for vcs information */
const char *lalInferenceVCSVersion = LALINFERENCE_VERSION;
const char *lalInferenceVCSId = LALINFERENCE_VCS_ID;
const char *lalInferenceVCSDate = LALINFERENCE_VCS_DATE;
const char *lalInferenceVCSBranch = LALINFERENCE_VCS_BRANCH;
const char *lalInferenceVCSTag = LALINFERENCE_VCS_TAG;
const char *lalInferenceVCSAuthor = LALINFERENCE_VCS_AUTHOR;
const char *lalInferenceVCSCommitter = LALINFERENCE_VCS_COMMITTER;
const char *lalInferenceVCSStatus = LALINFERENCE_VCS_STATUS;

/* global variables for vcs information - identable */
const char *lalInferenceVCSIdentId = LALINFERENCE_VCS_IDENT_ID;
const char *lalInferenceVCSIdentDate = LALINFERENCE_VCS_IDENT_DATE;
const char *lalInferenceVCSIdentBranch = LALINFERENCE_VCS_IDENT_BRANCH;
const char *lalInferenceVCSIdentTag = LALINFERENCE_VCS_IDENT_TAG;
const char *lalInferenceVCSIdentAuthor = LALINFERENCE_VCS_IDENT_AUTHOR;
const char *lalInferenceVCSIdentCommitter = LALINFERENCE_VCS_IDENT_COMMITTER;
const char *lalInferenceVCSIdentStatus = LALINFERENCE_VCS_IDENT_STATUS;

/* vcs information structure */
struct tagLALVCSInfo lalInferenceVCSInfo = { \
  LALINFERENCE_VCS_NAME, \
  LALINFERENCE_VERSION, \
  LALINFERENCE_VCS_ID, \
  LALINFERENCE_VCS_DATE, \
  LALINFERENCE_VCS_BRANCH, \
  LALINFERENCE_VCS_TAG, \
  LALINFERENCE_VCS_AUTHOR, \
  LALINFERENCE_VCS_COMMITTER, \
  LALINFERENCE_VCS_STATUS \
};

/*
 * vim: tw=0 ts=2 et
 */
