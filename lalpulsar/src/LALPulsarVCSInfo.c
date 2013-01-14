/*
 * LALPulsarVCSInfo.c - LALPulsar VCS Information
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

#include <lal/LALPulsarConfig.h>
#include <lal/LALPulsarVCSInfo.h>

/* global variables for vcs information */
const char *const lalPulsarVCSVersion = LALPULSAR_VERSION;
const char *const lalPulsarVCSId = LALPULSAR_VCS_ID;
const char *const lalPulsarVCSDate = LALPULSAR_VCS_DATE;
const char *const lalPulsarVCSBranch = LALPULSAR_VCS_BRANCH;
const char *const lalPulsarVCSTag = LALPULSAR_VCS_TAG;
const char *const lalPulsarVCSAuthor = LALPULSAR_VCS_AUTHOR;
const char *const lalPulsarVCSCommitter = LALPULSAR_VCS_COMMITTER;
const char *const lalPulsarVCSStatus = LALPULSAR_VCS_STATUS;

/* global variables for vcs information - identable */
const char *const lalPulsarVCSIdentId = LALPULSAR_VCS_IDENT_ID;
const char *const lalPulsarVCSIdentDate = LALPULSAR_VCS_IDENT_DATE;
const char *const lalPulsarVCSIdentBranch = LALPULSAR_VCS_IDENT_BRANCH;
const char *const lalPulsarVCSIdentTag = LALPULSAR_VCS_IDENT_TAG;
const char *const lalPulsarVCSIdentAuthor = LALPULSAR_VCS_IDENT_AUTHOR;
const char *const lalPulsarVCSIdentCommitter = LALPULSAR_VCS_IDENT_COMMITTER;
const char *const lalPulsarVCSIdentStatus = LALPULSAR_VCS_IDENT_STATUS;

/* vcs information structure */
const struct tagLALVCSInfo lalPulsarVCSInfo = { \
  LALPULSAR_VCS_NAME, \
  LALPULSAR_VERSION, \
  LALPULSAR_VCS_ID, \
  LALPULSAR_VCS_DATE, \
  LALPULSAR_VCS_BRANCH, \
  LALPULSAR_VCS_TAG, \
  LALPULSAR_VCS_AUTHOR, \
  LALPULSAR_VCS_COMMITTER, \
  LALPULSAR_VCS_STATUS \
};

/*
 * vim: tw=0 ts=2 et
 */
