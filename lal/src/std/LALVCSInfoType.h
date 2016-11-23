/*
 * Copyright (C) 2014, 2016 Karl Wette
 * Copyright (C) 2009-2013 Adam Mercer
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
 */

#ifndef _LALVCSINFOTYPE_H
#define _LALVCSINFOTYPE_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup LALVCSInfoType_h Header LALVCSInfoType.h
 * \ingroup lal_std
 * \author Adam Mercer, Karl Wette
 * \brief Contains routines for dealing with VCS and build information
 */
/*@{*/

/**
 * VCS and build information structure
 */
typedef struct tagLALVCSInfo
{
  const char *const name;		/**< Library name */
  const char *const version;		/**< Library version */
  const char *const vcsId;		/**< Identifier (e.g. git SHA1) of last commit */
  const char *const vcsDate;		/**< Committer date of last commit */
  const char *const vcsBranch;		/**< Branch of last commit */
  const char *const vcsTag;		/**< Tag of last commit */
  const char *const vcsAuthor;		/**< Author of last commit */
  const char *const vcsCommitter;	/**< Committer of last commit */
  const char *const vcsClean;		/**< (UN)CLEAN */
  const char *const vcsStatus;		/**< (UN)CLEAN: Status message */
  const char *const configureArgs;	/**< <tt>configure</tt> arguments */
  const char *const configureDate;	/**< <tt>configure</tt> date */
  const char *const buildDate;		/**< Build date */
} LALVCSInfo;

/**
 * <tt>NULL</tt>-terminated list of VCS and build information structures
 */
typedef const LALVCSInfo *const LALVCSInfoList[16];

/**
 * Compare two VCS information structures \p vcs1 and \p vcs2
 * \returns Zero if the structures are identical, non-zero otherwise
 */
int XLALVCSInfoCompare(const LALVCSInfo *vcs1, const LALVCSInfo *vcs2);

/**
 * Generate a multi-line string containing VCS and build information for a library and
 * its dependencies, as given in \p vcs_list. The verbosity of information contained in
 * the string is controlled by \p verbose. The string \p prefix is prepended to each line.
 */
char *XLALVCSInfoString(const LALVCSInfoList vcs_list, const int verbose, const char *prefix);

/*@}*/

#ifdef __cplusplus
}
#endif

#endif /* _LALVCSINFOTYPE_H */
