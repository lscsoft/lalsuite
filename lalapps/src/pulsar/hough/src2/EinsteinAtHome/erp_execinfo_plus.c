/***************************************************************************
 *   Copyright (C) 2009 by Oliver Bock                                     *
 *   oliver.bock[AT]aei.mpg.de                                             *
 *                                                                         *
 *   Derived from: addr2line.c (Ulrich Lauther, FSF Inc.)                  *
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

#include "sysdep.h"
#include "bfd.h"
#include "libiberty.h"
#include "demangle.h"

#include "erp_execinfo_plus.h"

#define MAX_NAME_LENGTH 512
#define MAX_ADDR_LENGTH 16

static bfd_boolean unwind_inlines; /* -i, unwind inlined functions. */
static bfd_boolean with_functions; /* -f, show function names.  */
static bfd_boolean do_demangle; /* -C, demangle names.  */
static bfd_boolean base_names; /* -s, strip directory names.  */

static asymbol **syms; /* Symbol table.  */
static FILE *ostream; /* output stream */

static int slurp_symtab(bfd *);
static void find_address_in_section(bfd *, asection *, void *);
static void translate_address(bfd *, const char *);

void bfd_nonfatal(const char *);
void bfd_fatal(const char *) ATTRIBUTE_NORETURN;
void report(const char *, va_list) ATTRIBUTE_PRINTF(1,0);
void fatal(const char *, ...) ATTRIBUTE_PRINTF_1 ATTRIBUTE_NORETURN;
void non_fatal(const char *, ...) ATTRIBUTE_PRINTF_1;
int set_default_bfd_target(void);
void list_matching_formats(char **);
off_t get_file_size(const char *);

/* Read in the symbol table.  */

static int slurp_symtab(bfd *abfd)
{
	long symcount;
	unsigned int size;

	if ((bfd_get_file_flags (abfd)& HAS_SYMS) == 0)
		return 0;

	symcount = bfd_read_minisymbols (abfd, FALSE, (void **) &syms, &size);
	if (symcount == 0)
		symcount = bfd_read_minisymbols (abfd, TRUE /* dynamic */, (void **) &syms, &size);

	if (symcount < 0) {
		bfd_fatal(bfd_get_filename (abfd));
		return -1;
	}
	return 0;
}

/* These global variables are used to pass information between
 translate_addresses and find_address_in_section.  */

static bfd_vma pc;
static const char *filename;
static const char *functionname;
static unsigned int line;
static bfd_boolean found;

/* Look for an address in a section.  This is called via
 bfd_map_over_sections.  */

static void find_address_in_section(bfd *abfd, asection *section, void *data ATTRIBUTE_UNUSED)
{
	bfd_vma vma;
	bfd_size_type size;

	if (found)
		return;

	if ((bfd_get_section_flags (abfd, section)& SEC_ALLOC) == 0)
		return;

	vma = bfd_get_section_vma (abfd, section);
	if (pc < vma)
		return;

	size = bfd_get_section_size (section);
	if (pc >= vma + size)
		return;

	found = bfd_find_nearest_line (abfd, section, syms, pc - vma,
			&filename, &functionname, &line);
}

/* Translate hexadecimal address into
 file_name:line_number and optionally function name.  */

static void translate_address(bfd *abfd, const char *address)
{
	pc = bfd_scan_vma(address, NULL, 16);

	found = FALSE;

	bfd_map_over_sections(abfd, find_address_in_section, NULL);

	if (found) {
		do {
			const char *fnName;
			if (with_functions) {
				char *alloc = NULL;

				fnName = functionname;
				if (fnName == NULL || *fnName == '\0')
					fnName = "unknown";
				else if (do_demangle) {
					alloc = bfd_demangle(abfd, fnName, DMGL_ANSI | DMGL_PARAMS);
					if (alloc != NULL)
						fnName = alloc;
				}

				if (alloc != NULL)
					free(alloc);
			} else {
				fnName = "disabled";
			}

			if (base_names && filename != NULL) {
				char *h;

				h = strrchr(filename, '/');
				if (h != NULL)
					filename = h + 1;
			}

			fprintf(ostream, "\tSource file: %s (Function: %s / Line: %u)\n",
					filename ? filename : "unknown", fnName, line);
			if (!unwind_inlines)
				found = FALSE;
			else
				found = bfd_find_inliner_info (abfd, &filename, &functionname, &line);
		} while (found);

	}

	fflush(ostream);
}

/* Process a file */

static int process_file(const char *file_name, const char *address)
{
	bfd *abfd;
	char **matching;

	if (get_file_size(file_name) < 1)
		return -1;

	abfd = bfd_openr(file_name, NULL);
	if (abfd == NULL) {
		bfd_fatal(file_name);
		return -1;
	}

	if (bfd_check_format(abfd, bfd_archive)) {
		fatal(_("%s: cannot get addresses from archive"), file_name);
		return -1;
	}

	if (!bfd_check_format_matches(abfd, bfd_object, &matching)) {
		bfd_nonfatal(bfd_get_filename (abfd));
		if (bfd_get_error() == bfd_error_file_ambiguously_recognized) {
			list_matching_formats(matching);
			free(matching);
		}
		return 1;
	}

	if (slurp_symtab(abfd) != 0) {
		return -1;
	}

	translate_address(abfd, address);

	if (syms != NULL) {
		free(syms);
		syms = NULL;
	}

	bfd_close(abfd);

	return 0;
}

int backtrace_symbols_fd_plus(const char *const *symbols, int size, int fd)
{
	int i;

	/* sanity check */
	if (fd <= 0) {
		fd = fileno(stderr);
	}

	/* get file handle */
	ostream = fdopen(fd, "a");
	if (!ostream) {
		fprintf(stderr, "Error: couldn't open output stream!");
		return -1;
	}

	/* sanity check */
	if (!symbols || size <= 0) {
		fprintf(ostream, "Error: invalid or empty symbol list encountered!");
		return -1;
	}

	/* default config */
	base_names = TRUE;
	do_demangle = TRUE;
	unwind_inlines = TRUE;
	with_functions = TRUE;

	/* init libbfd */
	int retval = 0;
	bfd_init();
	retval = set_default_bfd_target();

	if (retval != 0) {
		return retval;
	}

	for (i = 0; i < size; ++i) {
		/* content buffers */
		char file_name[MAX_NAME_LENGTH +1] = { '\0' };
		char func_offset[MAX_NAME_LENGTH +1] = { '\0' };
		char address[MAX_NAME_LENGTH +1] = { '\0' };

		/* find file name delimiter */
		char *fnEnd = strrchr(symbols[i], ' ');
		/* find funtion+offset delimiters */
		char *firstP = strrchr(symbols[i], '(');
		char *secondP = strrchr(symbols[i], ')');
		/* find address delimiters */
		char *firstB = strrchr(symbols[i], '[');
		char *secondB = strrchr(symbols[i], ']');

		/* get file name */
		size_t length = 0;
		if (firstP != NULL) {
			length = firstP - symbols[i];
		} else if (fnEnd != NULL) {
			length = fnEnd - symbols[i];
		} else {
			fprintf(
					ostream,
					"Error: couldn't parse symbol information for file name: %s",
					symbols[i]);
			retval = 1;
			break;
		}
		if (length > MAX_NAME_LENGTH) {
			length = MAX_NAME_LENGTH;
		}
		strncpy(file_name, symbols[i], length);

		/* get address */
		if (firstB != NULL && secondB != NULL) {
			length = secondB - firstB;
			if (length > MAX_ADDR_LENGTH) {
				length = MAX_ADDR_LENGTH;
			}
			strncpy(address, firstB + 1, length - 1);
		} else {
			fprintf(ostream,
					"Error: couldn't parse symbol information for address: %s",
					symbols[i]);
			retval = 1;
			break;
		}

		/* get function-offset if available */
		if (firstP != NULL && secondP != NULL) {
			length = secondP - firstP;
			if (length > MAX_NAME_LENGTH) {
				length = MAX_NAME_LENGTH;
			}
			strncpy(func_offset, firstP + 1, length - 1);
		}

		/* finally, get symbol information */
		if (strlen(file_name) > 0 && strlen(file_name) > 0) {
			fprintf(ostream, "Frame %d:\n\tBinary file: %s (%s)\n", size - i, file_name, address);
			if (strlen(func_offset) > 0) {
				fprintf(ostream, "\tOffset info: %s\n", func_offset);
			}
			retval = process_file(file_name, address);
		}
	}

	return retval;
}

void bfd_nonfatal(const char *string)
{
	const char *errmsg = bfd_errmsg(bfd_get_error());

	if (string)
		fprintf(ostream, "%s: %s\n", string, errmsg);
	else
		fprintf(ostream, "%s\n", errmsg);
}

void bfd_fatal(const char *string)
{
	bfd_nonfatal(string);
}

void report(const char * format, va_list args)
{
	vfprintf(ostream, format, args);
	putc ('\n', ostream);
}

void fatal VPARAMS ((const char *format, ...))
{
	VA_OPEN (args, format);
	VA_FIXEDARG (args, const char *, format);

	report(format, args);
	VA_CLOSE (args);
}

void non_fatal VPARAMS ((const char *format, ...))
{
	VA_OPEN (args, format);
	VA_FIXEDARG (args, const char *, format);

	report(format, args);
	VA_CLOSE (args);
}

/* Returns the size of the named file.  If the file does not
 exist, or if it is not a real file, then a suitable non-fatal
 error message is printed and zero is returned.  */

off_t get_file_size(const char * file_name)
{
	struct stat statbuf;

	if (stat(file_name, &statbuf) < 0) {
		if (errno == ENOENT)
			non_fatal(_("'%s': No such file"), file_name);
		else
			non_fatal(_("Warning: could not locate '%s'.  reason: %s"),
					file_name, strerror(errno));
	} else if (!S_ISREG(statbuf.st_mode))
		non_fatal(_("Warning: '%s' is not an ordinary file"), file_name);
	else
		return statbuf.st_size;

	return 0;
}

/* After a FALSE return from bfd_check_format_matches with
 bfd_get_error () == bfd_error_file_ambiguously_recognized, print
 the possible matching targets.  */

void list_matching_formats(char **p)
{
	fprintf(ostream, _("Matching formats:"));
	while (*p)
		fprintf(ostream, " %s", *p++);
	fputc('\n', ostream);
}

/* Set the default BFD target based on the configured target.  Doing
 this permits the binutils to be configured for a particular target,
 and linked against a shared BFD library which was configured for a
 different target.  */

int set_default_bfd_target(void)
{
	/* The macro TARGET is defined by Makefile.  */
	const char *target = TARGET;

	if (!bfd_set_default_target(target)) {
		fatal(_("can't set BFD default target to `%s': %s"), target,
				bfd_errmsg(bfd_get_error()));
		return -1;
	}
}
