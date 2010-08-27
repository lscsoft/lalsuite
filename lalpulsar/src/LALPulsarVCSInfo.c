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
const char *lalPulsarVCSVersion = LALPULSAR_VERSION;
const char *lalPulsarVCSId = LALPULSAR_VCS_ID;
const char *lalPulsarVCSDate = LALPULSAR_VCS_DATE;
const char *lalPulsarVCSBranch = LALPULSAR_VCS_BRANCH;
const char *lalPulsarVCSTag = LALPULSAR_VCS_TAG;
const char *lalPulsarVCSAuthor = LALPULSAR_VCS_AUTHOR;
const char *lalPulsarVCSCommitter = LALPULSAR_VCS_COMMITTER;
const char *lalPulsarVCSStatus = LALPULSAR_VCS_STATUS;

/* global variables for vcs information - identable */
const char *lalPulsarVCSIdentId = LALPULSAR_VCS_IDENT_ID;
const char *lalPulsarVCSIdentDate = LALPULSAR_VCS_IDENT_DATE;
const char *lalPulsarVCSIdentBranch = LALPULSAR_VCS_IDENT_BRANCH;
const char *lalPulsarVCSIdentTag = LALPULSAR_VCS_IDENT_TAG;
const char *lalPulsarVCSIdentAuthor = LALPULSAR_VCS_IDENT_AUTHOR;
const char *lalPulsarVCSIdentCommitter = LALPULSAR_VCS_IDENT_COMMITTER;
const char *lalPulsarVCSIdentStatus = LALPULSAR_VCS_IDENT_STATUS;

/* vcs information structure */
struct tagLALVCSInfo lalPulsarVCSInfo = { \
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
