/*
 * LALFrameVCSInfo.c - LALFrame VCS Information
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

#include <lal/LALFrameConfig.h>
#include <lal/LALFrameVCSInfo.h>

/* global variables for vcs information */
const char *const lalFrameVCSVersion = LALFRAME_VERSION;
const char *const lalFrameVCSId = LALFRAME_VCS_ID;
const char *const lalFrameVCSDate = LALFRAME_VCS_DATE;
const char *const lalFrameVCSBranch = LALFRAME_VCS_BRANCH;
const char *const lalFrameVCSTag = LALFRAME_VCS_TAG;
const char *const lalFrameVCSAuthor = LALFRAME_VCS_AUTHOR;
const char *const lalFrameVCSCommitter = LALFRAME_VCS_COMMITTER;
const char *const lalFrameVCSStatus = LALFRAME_VCS_STATUS;

/* global variables for vcs information - identable */
const char *const lalFrameVCSIdentId = LALFRAME_VCS_IDENT_ID;
const char *const lalFrameVCSIdentDate = LALFRAME_VCS_IDENT_DATE;
const char *const lalFrameVCSIdentBranch = LALFRAME_VCS_IDENT_BRANCH;
const char *const lalFrameVCSIdentTag = LALFRAME_VCS_IDENT_TAG;
const char *const lalFrameVCSIdentAuthor = LALFRAME_VCS_IDENT_AUTHOR;
const char *const lalFrameVCSIdentCommitter = LALFRAME_VCS_IDENT_COMMITTER;
const char *const lalFrameVCSIdentStatus = LALFRAME_VCS_IDENT_STATUS;

/* vcs information structure */
const struct tagLALVCSInfo lalFrameVCSInfo = { LALFRAME_VCS_NAME,
    LALFRAME_VERSION,
    LALFRAME_VCS_ID,
    LALFRAME_VCS_DATE,
    LALFRAME_VCS_BRANCH,
    LALFRAME_VCS_TAG,
    LALFRAME_VCS_AUTHOR,
    LALFRAME_VCS_COMMITTER,
    LALFRAME_VCS_STATUS
};

/*
 * vim: tw=0 ts=2 et
 */
