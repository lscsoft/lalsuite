/*
*  Copyright (C) 2013 Jolien Creighton
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <limits.h>
#include <string.h>
#include <unistd.h>
#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/LALSimReadData.h>

#ifndef PAGESIZE
#ifdef _SC_PAGE_SIZE
#define PAGESIZE _SC_PAGE_SIZE
#else
#define PAGESIZE 1024
#endif
#endif
#ifndef LINE_MAX
#ifdef _SC_LINE_MAX
#define LINE_MAX _SC_LINE_MAX
#else
#define LINE_MAX 1024
#endif
#endif


/**
 * @brief Opens a specified data file, searching default path if necessary.
 * @details Opens a data file for input with a specified path name.
 * If the path name is an absolute path then this specific file is opened;
 * otherwise, search for the file in paths given in the environment variable
 * LALSIM_DATA_PATH, and finally search in the installed PKG_DATA_DIR path.
 * @param[in] fname The path of the file to open.
 * @return A pointer to a LALFILE structure or NULL on failure.
 */
LALFILE *XLALSimReadDataFileOpen(const char *fname)
{
    const char *pkgdatadir = PKG_DATA_DIR;
    char path[PATH_MAX] = "";
    LALFILE *fp;

    if (strchr(fname, '/')) {   /* a specific path is given */
        if (realpath(fname, path) == NULL)
            XLAL_ERROR_NULL(XLAL_EIO, "Unresolvable path %s\n", path);
    } else {
        /* unspecific path given: use LALSIM_DATA_PATH environment */
        char *env = getenv("LALSIM_DATA_PATH");
        char *str;
        char *dir;
        env = str = XLALStringDuplicate(env ? env : ":");
        while ((dir = strsep(&str, ":"))) {
            if (strlen(dir))
                snprintf(path, sizeof(path), "%s/%s", dir, fname);
            else        /* use default path */
                snprintf(path, sizeof(path), "%s/%s", pkgdatadir, fname);
            if (access(path, R_OK) == 0)        /* found it! */
                break;
            *path = 0;
        }
        XLALFree(env);
    }
    if (!*path) /* could not find file */
        XLAL_ERROR_NULL(XLAL_EIO, "Could not find data file %s\n", fname);
    fp = XLALFileOpenRead(path);
    if (!fp)    /* open failure */
        XLAL_ERROR_NULL(XLAL_EIO, "Could not open data file %s\n", path);
    return fp;
}


/**
 * @brief Read a two-column data file.
 * @details Read a data file containing two whitespace separated columns
 * of data and create two arrays containing the data in each column.
 * If any line begins with the character '#' then it is ignored.
 * @param[out] xdat The x-data stored in the first column.
 * @param[out] ydat The y-data stored in the second column.
 * @param fp Pointer to a LALFILE structure opened for input.
 * @return The number of data points read or <0 if an error occurs.
 */
size_t XLALSimReadDataFile2Col(double **xdat, double **ydat, LALFILE * fp)
{
    char line[LINE_MAX];
    size_t size = PAGESIZE;
    size_t lnum = 0;
    size_t npts;
    *xdat = XLALMalloc(size * sizeof(**xdat));
    *ydat = XLALMalloc(size * sizeof(**ydat));
    npts = 0;
    while (XLALFileGets(line, sizeof(line), fp)) {
        ++lnum;
        if (strchr(line, '\n') == NULL) {   /* line too long */
            XLALFree(*xdat);
            XLALFree(*ydat);
            XLAL_ERROR(XLAL_EIO, "Line %zd too long\n", lnum);
        }
        if (*line == '#')       /* ignore lines beginning with a '#' */
            continue;
        if (sscanf(line, "%lf %lf", *xdat + npts, *ydat + npts) != 2) {
            XLALFree(*xdat);
            XLALFree(*ydat);
            XLAL_ERROR(XLAL_EIO, "Line %zd malformed\n", lnum);
        }
        if (++npts == size) {
            size += PAGESIZE;
            *xdat = XLALRealloc(*xdat, size * sizeof(**xdat));
            *ydat = XLALRealloc(*ydat, size * sizeof(**ydat));
        }
    }
    *xdat = XLALRealloc(*xdat, npts * sizeof(**xdat));
    *ydat = XLALRealloc(*ydat, npts * sizeof(**ydat));
    return npts;
}


/**
 * @brief Read a multi-column data file.
 * @details Read a data file containing multiple whitespace separated columns
 * of data and create an array containing the data.
 * If any line begins with the character '#' then it is ignored.
 * The data is stored in the array in row-major format so that the data
 * sample on row @c i (beginning with zero) and column @c j (beginning with
 * zero) is found as the element <tt>[i * ncol + j]</tt> where @c ncol is the
 * number of columns.
 * @param[out] data The data stored in row-major order.
 * @param[out] ncol The number of columns in the data file.
 * @param fp Pointer to a LALFILE structure opened for input.
 * @return The number of rows read or (size_t)(-1) if an error occurs.
 */
size_t XLALSimReadDataFileNCol(double **data, size_t *ncol, LALFILE *fp)
{
    char line[LINE_MAX];
    size_t page = PAGESIZE;
    size_t size = 0;
    size_t lnum = 0;
    size_t nrow = 0;

    *data = NULL;
    *ncol = 0;
    while (XLALFileGets(line, sizeof(line), fp)) {
        char *s;
        char *endp;
        size_t col;

        ++lnum;

        if (strchr(line, '\n') == NULL) {   /* line too long */
            XLALFree(*data);
            XLAL_ERROR(XLAL_EIO, "Line %zd too long\n", lnum);
        }

        if (*line == '#')       /* ignore lines beginning with '#' */
            continue;

        if (*ncol == 0) {       /* count columns on first line */
            endp = line;
            while (1) {
                s = endp;
                /* work around bug in glibc < 2.16
                 * http://sourceware.org/bugzilla/show_bug.cgi?id=13970 */
                double v = strtod(s, &endp);
                (void)v;
                if (s == endp || *endp == '\0')
                    break;
                ++*ncol;
            }
            if (*ncol == 0) {
                XLALFree(*data);
                XLAL_ERROR(XLAL_EIO, "Line %zd malformed\n", lnum);
            }
        }

        if (nrow == size) {     /* allocate more memory for data */
            size += page;
            *data = XLALRealloc(*data, *ncol * page * sizeof(**data));
        }

        /* scan line for data values in each column */
        endp = line;
        for (col = 0; col < *ncol; ++col) {
            s = endp;
            (*data)[*ncol * nrow + col] = strtod(s, &endp);
            if (s == endp || *endp == '\0') {
                XLALFree(*data);
                XLAL_ERROR(XLAL_EIO, "Line %zd malformed\n", lnum);
            }
        }
        
        ++nrow;
    }

    *data = XLALRealloc(*data, *ncol * nrow * sizeof(**data));

    return nrow;
}
