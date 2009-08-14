/***************************************************************************
 *   Copyright (C) 2009 by Oliver Bock                                     *
 *   oliver.bock[AT]aei.mpg.de                                             *
 *                                                                         *
 *   This file is part of Einstein@Home (Radio Pulsar Edition).            *
 *                                                                         *
 *   Einstein@Home is free software: you can redistribute it and/or modify *
 *   it under the terms of the GNU General Public License as published     *
 *   by the Free Software Foundation, version 2 of the License.            *
 *                                                                         *
 *   Einstein@Home is distributed in the hope that it will be useful,      *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with Einstein@Home. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                         *
 ***************************************************************************/

#ifndef ERP_EXECINFO_PLUS_H
#define ERP_EXECINFO_PLUS_H

#ifdef __cplusplus
extern "C" {
#endif

	/**
	 * Takes the string list of symbols returned by glibc's
	 * backtrace_symbols() and prints its content with additional
	 * information to the given file descriptor.
	 *
	 * @return: 0 if success, <0 if error, >0 if warning
	 */
	int backtrace_symbols_fd_plus(const char *const *symbols, int size, int fd);

#ifdef __cplusplus
}
#endif

#endif
