/*
*  Copyright (C) 2007 Jolien Creighton
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <config.h>

#ifdef HAVE_REGEX_H
/* hack to get around problem in regex.h */
#ifdef __GLIBC__
#ifdef __restrict_arr
#undef __restrict_arr
#endif
#define __restrict_arr
#endif
#include <regex.h>
#endif

#ifdef HAVE_GLOB_H
#include <glob.h>
#endif

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/Sort.h>
#include <lal/LALCache.h>
#include <lal/FileIO.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

static int XLALCacheFileReadRow(char *s, size_t len, LALFILE * fp,
                                int *line)
{
    while (1) {
        if (!XLALFileGets(s, len, fp))
            return 0;   /* finished code */
        ++(*line);
        if (!strchr(s, '\n')) {
            if (XLALFileEOF(fp)) {      /* missing final newline */
                XLAL_PRINT_WARNING("Missing newline on line %d", *line);
                return 0;
            }
            /* line is too long */
            XLAL_ERROR(XLAL_EIO, "Line %d too long", *line);
        }
        if (*s != '#')
            break;
    }
    return 1;   /* continue code */
}

/* counts the rows of the file */
static int XLALCacheFileCountRows(LALFILE * fp)
{
    char s[2 * FILENAME_MAX];
    int line = 0;
    int n = 0;
    int c;
    while ((c = XLALCacheFileReadRow(s, sizeof(s), fp, &line)))
        if (c < 0)
            return -1;
        else
            ++n;
    return n;
}

/* kind of like strtok -- modifies buffer */
static char *XLALCacheFileNextField(char **ps)
{
    const char *sepstr = " \t\n";
    char *s = *ps;
    char *field;
    size_t len;
    len = strcspn(s, sepstr);
    if (!len)
        return NULL;
    field = s;
    s += len;
    *s++ = 0;
    s += strspn(s, sepstr);
    *ps = s;
    return field;
}


/* Parses a row into appropriate fields in the cache entry */
/* NOTE: perfectly happy if there is too many rows! */
static int XLALCacheFileParseEntry(struct tagLALCacheEntry *entry, char *s)
{
    char *f1, *f2, *f3, *f4, *f5;
    if ((f1 = XLALCacheFileNextField(&s))
        && (f2 = XLALCacheFileNextField(&s))
        && (f3 = XLALCacheFileNextField(&s))
        && (f4 = XLALCacheFileNextField(&s))
        && (f5 = XLALCacheFileNextField(&s))) {
        entry->src = strcmp(f1, "-") ? XLALStringDuplicate(f1) : NULL;
        entry->dsc = strcmp(f2, "-") ? XLALStringDuplicate(f2) : NULL;
        entry->url = strcmp(f5, "-") ? XLALStringDuplicate(f5) : NULL;
        if (strcmp(f3, "-")) {
            if (strspn(f3, "0123456789") != strlen(f3))
                XLAL_ERROR(XLAL_EIO, "Invalid content in field 3 \"%s\"",
                           f3);
            entry->t0 = atoi(f3);
        } else
            entry->t0 = 0;
        if (strcmp(f4, "-")) {
            if (strspn(f4, "0123456789") != strlen(f4))
                XLAL_ERROR(XLAL_EIO, "Invalid content in field 4 \"%s\"",
                           f4);
            entry->dt = atoi(f4);
        } else
            entry->dt = 0;
        return 0;
    }
    XLAL_ERROR(XLAL_EIO, "Wrong number of fields");     /* wrong field count */
}

static int XLALCacheEntryCopy(LALCacheEntry * dst,
                              const LALCacheEntry * src)
{
    if (!src || !dst)
        XLAL_ERROR(XLAL_EFAULT);
    dst->src = XLALStringDuplicate(src->src);
    dst->dsc = XLALStringDuplicate(src->dsc);
    dst->url = XLALStringDuplicate(src->url);
    dst->t0 = src->t0;
    dst->dt = src->dt;
    return 0;
}

LALCache *XLALCreateCache(UINT4 length)
{
    LALCache *cache;
    cache = XLALCalloc(1, sizeof(*cache));
    if (!cache)
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    if (length) {
        cache->list = XLALCalloc(length, sizeof(*cache->list));
        if (!cache->list) {
            XLALFree(cache);
            XLAL_ERROR_NULL(XLAL_ENOMEM);
        }
    } else
        XLAL_PRINT_WARNING("Creating a zero-length cache");
    cache->length = length;
    return cache;
}

void XLALDestroyCache(LALCache * cache)
{
    if (cache) {
        UINT4 i;
        for (i = 0; i < cache->length; ++i) {
            XLALFree(cache->list[i].src);
            XLALFree(cache->list[i].dsc);
            XLALFree(cache->list[i].url);
        }
        XLALFree(cache->list);
        XLALFree(cache);
    }
    return;
}

LALCache *XLALCacheDuplicate(const LALCache * cache)
{
    LALCache *duplicate = NULL;
    if (cache) {
        UINT4 i;
        duplicate = XLALCreateCache(cache->length);
        if (!duplicate)
            XLAL_ERROR_NULL(XLAL_EFUNC);
        for (i = 0; i < cache->length; ++i)
            XLALCacheEntryCopy(duplicate->list + i, cache->list + i);
    }
    return duplicate;
}

LALCache *XLALCacheMerge(const LALCache * cache1, const LALCache * cache2)
{
    LALCache *cache = NULL;
    LALCacheEntry *entry;
    UINT4 length;
    UINT4 i;
    if (!cache2)
        return XLALCacheDuplicate(cache1);
    if (!cache1)
        return XLALCacheDuplicate(cache2);
    length = cache1->length + cache2->length;
    cache = XLALCreateCache(length);
    if (!cache)
        XLAL_ERROR_NULL(XLAL_EFUNC);
    entry = cache->list;
    for (i = 0; i < cache1->length; ++i)
        XLALCacheEntryCopy(entry++, cache1->list + i);
    for (i = 0; i < cache2->length; ++i)
        XLALCacheEntryCopy(entry++, cache2->list + i);

    XLALCacheSort(cache);
    return cache;
}

LALCache *XLALCacheFileRead(LALFILE * fp)
{
    LALCache *cache;
    char s[2 * FILENAME_MAX];
    int line = 0;
    int n;
    int i;
    if (!fp)
        XLAL_ERROR_NULL(XLAL_EFAULT);
    n = XLALCacheFileCountRows(fp);
    if (n < 0)
        XLAL_ERROR_NULL(XLAL_EFUNC);
    XLALFileRewind(fp);
    cache = XLALCreateCache(n);
    if (!cache)
        XLAL_ERROR_NULL(XLAL_EFUNC);
    for (i = 0; i < n; ++i)
        if (XLALCacheFileReadRow(s, sizeof(s), fp, &line) != 1
            || XLALCacheFileParseEntry(&cache->list[i], s) < 0) {
            XLALDestroyCache(cache);
            XLAL_ERROR_NULL(XLAL_EFUNC, "Error reading row %s on line %s",
                            i + 1, line);
        }
    XLALCacheSort(cache);
    return cache;
}

LALCache *XLALCacheImport(const char *fname)
{
    LALCache *cache;
    LALFILE *fp;
    fp = XLALFileOpenRead(fname);
    if (!fp)
        XLAL_ERROR_NULL(XLAL_EIO);
    cache = XLALCacheFileRead(fp);
    XLALFileClose(fp);
    if (!fp)
        XLAL_ERROR_NULL(XLAL_EFUNC);
    return cache;
}

#ifdef HAVE_GLOB_H      /* only use this if globbing is supported */
static int XLALCacheFilenameParseEntry(LALCacheEntry * entry,
                                       const char *path)
{
    char src[FILENAME_MAX];
    char dsc[FILENAME_MAX];
    const char *base;
    INT4 t0;
    INT4 dt;
    int c;

    /* basename */
    base = strrchr(path, '/');
    base = base ? base + 1 : path;

    /* construct url */
    entry->url = XLALStringDuplicate("file://localhost");
    if (!entry->url)
        XLAL_ERROR(XLAL_EFUNC);
    if (*path == '/') { /* absolute path */
        entry->url = XLALStringAppend(entry->url, path);
    } else {    /* relative path */
        char cwd[FILENAME_MAX];
#if 0
        char *ptr;
        ptr = getcwd(cwd, sizeof(cwd) - 1);
#endif
        XLALStringConcatenate(cwd, "/", sizeof(cwd));
        XLALStringConcatenate(cwd, path, sizeof(cwd));
        entry->url = XLALStringAppend(entry->url, cwd);
    }
    if (!entry->url)
        XLAL_ERROR(XLAL_EFUNC);

    /* extract src, dsc, t0, and dt from file name */
    c = sscanf(base, "%[a-zA-Z0-9_+#]-%[a-zA-Z0-9_+#]-%d-%d", src, dsc,
               &t0, &dt);
    if (c == 4) {       /* expected format */
        entry->src = XLALStringDuplicate(src);
        entry->dsc = XLALStringDuplicate(dsc);
        entry->t0 = t0;
        entry->dt = dt;
        if (!entry->src || !entry->dsc) {
            XLALFree(entry->src);
            XLALFree(entry->dsc);
            XLALFree(entry->url);
            XLAL_ERROR(XLAL_EFUNC);
        }
    }
    return 0;
}

#endif /* HAVE_GLOB_H */

LALCache *XLALCacheGlob(const char *dirstr, const char *fnptrn)
{
#ifdef HAVE_GLOB_H
    LALCache *cache;
    int globflags = 0;
    glob_t g;
    size_t i;

    fnptrn = fnptrn ? fnptrn : "*";
    dirstr = dirstr ? dirstr : ".";

    if (fnptrn[0]
        && (fnptrn[0] == '/' || (fnptrn[0] == '.' && fnptrn[1]
                                 && (fnptrn[1] == '/'
                                     || (fnptrn[1] == '.'
                                         && fnptrn[2] == '/')))))
        glob(fnptrn, globflags, NULL, &g);
    else {      /* prepend path from dirname */
        char path[FILENAME_MAX];
        char dirname[FILENAME_MAX];
        char *nextdir;
        strncpy(dirname, dirstr, sizeof(dirname) - 1);
        dirname[sizeof(dirname) - 1] = 0;
        do {
            if ((nextdir = strchr(dirname, ':')))
                *nextdir++ = 0;
            snprintf(path, sizeof(path), "%s/%s", *dirname ? dirname : ".",
                     fnptrn);
            glob(path, globflags, NULL, &g);
            fnptrn = path;
            globflags |= GLOB_APPEND;
        } while (nextdir);
    }

    if (!g.gl_pathc)
        XLAL_ERROR_NULL(XLAL_EIO, "No matching files found in %s\n",
                        fnptrn);
    cache = XLALCreateCache(g.gl_pathc);
    if (!cache) {
        globfree(&g);
        XLAL_ERROR_NULL(XLAL_EFUNC);
    }
    /* copy file names */
    for (i = 0; i < g.gl_pathc; ++i) {
        LALCacheEntry *entry = cache->list + i;
        if (0 > XLALCacheFilenameParseEntry(entry, g.gl_pathv[i])) {
            globfree(&g);
            XLALDestroyCache(cache);
            XLAL_ERROR_NULL(XLAL_EFUNC);
        }
    }

    globfree(&g);
    XLALCacheSort(cache);
    return cache;
#else /* no globbing: unsupported */
    fnptrn = NULL;
    dirstr = NULL;
    XLAL_ERROR_NULL(XLAL_EFAILED,
                    "Glob is unsupported on non-posix system");
#endif
}

int XLALCacheFileWrite(LALFILE * fp, const LALCache * cache)
{
    UINT4 i;
    if (!fp || !cache)
        XLAL_ERROR(XLAL_EFAULT);
    for (i = 0; i < cache->length; ++i) {
        XLALFilePrintf(fp, "%s\t",
                       cache->list[i].src ? cache->list[i].src : "-");
        XLALFilePrintf(fp, "%s\t",
                       cache->list[i].dsc ? cache->list[i].dsc : "-");
        if (cache->list[i].t0 > 0)
            XLALFilePrintf(fp, "%d\t", cache->list[i].t0);
        else
            XLALFilePrintf(fp, "-\t");
        if (cache->list[i].dt > 0)
            XLALFilePrintf(fp, "%d\t", cache->list[i].dt);
        else
            XLALFilePrintf(fp, "-\t");
        XLALFilePrintf(fp, "%s\n",
                       cache->list[i].url ? cache->list[i].url : "-");
    }
    return 0;
}


static int XLALCacheCompareSource(void UNUSED * p, const void *p1,
                                  const void *p2)
{
    const char *s1 = ((const struct tagLALCacheEntry *) p1)->src;
    const char *s2 = ((const struct tagLALCacheEntry *) p2)->src;
    return strcmp(s1 ? s1 : "", s2 ? s2 : "");
}

static int XLALCacheCompareDescription(void UNUSED * p, const void *p1,
                                       const void *p2)
{
    const char *s1 = ((const struct tagLALCacheEntry *) p1)->dsc;
    const char *s2 = ((const struct tagLALCacheEntry *) p2)->dsc;
    return strcmp(s1 ? s1 : "", s2 ? s2 : "");
}

static int XLALCacheCompareStartTime(void UNUSED * p, const void *p1,
                                     const void *p2)
{
    double t1 = ((const struct tagLALCacheEntry *) p1)->t0;
    double t2 = ((const struct tagLALCacheEntry *) p2)->t0;
    return (t1 > t2) - (t1 < t2);
}

static int XLALCacheCompareDuration(void UNUSED * p, const void *p1,
                                    const void *p2)
{
    double t1 = ((const struct tagLALCacheEntry *) p1)->dt;
    double t2 = ((const struct tagLALCacheEntry *) p2)->dt;
    return (t1 > t2) - (t1 < t2);
}

static int XLALCacheCompareURL(void UNUSED * p, const void *p1,
                               const void *p2)
{
    const char *s1 = ((const struct tagLALCacheEntry *) p1)->url;
    const char *s2 = ((const struct tagLALCacheEntry *) p2)->url;
    return strcmp(s1 ? s1 : "", s2 ? s2 : "");
}

static int XLALCacheCompareEntries(void *p, const void *p1, const void *p2)
{
    int c;
    if ((c = XLALCacheCompareSource(p, p1, p2)))
        return c;
    if ((c = XLALCacheCompareDescription(p, p1, p2)))
        return c;
    if ((c = XLALCacheCompareStartTime(p, p1, p2)))
        return c;
    if ((c = XLALCacheCompareDuration(p, p1, p2)))
        return c;
    return XLALCacheCompareURL(p, p1, p2);
}

static int XLALCacheCompareEntryMetadata(void *p, const void *p1,
                                         const void *p2)
{
    int c;
    if ((c = XLALCacheCompareSource(p, p1, p2)))
        return c;
    if ((c = XLALCacheCompareDescription(p, p1, p2)))
        return c;
    if ((c = XLALCacheCompareStartTime(p, p1, p2)))
        return c;
    if ((c = XLALCacheCompareDuration(p, p1, p2)))
        return c;
    return 0;
}


int XLALCacheSort(LALCache * cache)
{
    if (!cache)
        XLAL_ERROR(XLAL_EFAULT);
    return XLALHeapSort(cache->list, cache->length, sizeof(*cache->list),
                        NULL, XLALCacheCompareEntries);
}

int XLALCacheUniq(LALCache * cache)
{
    LALCacheEntry *list;
    UINT4 i, j;

    if (!cache)
        XLAL_ERROR(XLAL_EFAULT);
    list = cache->list;
    i = 0;
    for (j = 0; j < cache->length; ++j) {
        LALCacheEntry swap;
        swap = list[i];
        list[i] = list[j];
        list[j] = swap;
        if (j + 1 == cache->length
            || XLALCacheCompareEntryMetadata(NULL, list + i, list + j + 1))
            ++i;
    }
    for (j = i; j < cache->length; ++j) {
        XLALFree(list[j].src);
        XLALFree(list[j].dsc);
        XLALFree(list[j].url);
    }
    if (i != cache->length) {
        cache->length = i;
        cache->list = XLALRealloc(cache->list, i * sizeof(*cache->list));
        if (!cache->list)
            XLAL_ERROR(XLAL_EFUNC);
    }
    return 0;
}

#ifdef HAVE_REGEX_H
static int XLALCacheEntryMatch(const LALCacheEntry * entry, INT4 t0,
                               INT4 t1, regex_t * srcreg, regex_t * dscreg,
                               regex_t * urlreg)
{
    if (t1 > 0 && !(entry->t0 < t1))
        return 0;
    if (t0 > 0 && !(entry->t0 + entry->dt > t0))
        return 0;
    if (srcreg)
        if (!entry->src || regexec(srcreg, entry->src, 0, NULL, 0))
            return 0;
    if (dscreg)
        if (!entry->dsc || regexec(dscreg, entry->dsc, 0, NULL, 0))
            return 0;
    if (urlreg)
        if (!entry->url || regexec(urlreg, entry->url, 0, NULL, 0))
            return 0;
    return 1;   /* matches */
}

#else /* HAVE_REGEX_H undefined */
/* can only attempt to match time range */
static int XLALCacheEntryMatchTime(const LALCacheEntry * entry, INT4 t0,
                                   INT4 t1)
{
    if (t1 > 0 && !(entry->t0 < t1))
        return 0;
    if (t0 > 0 && !(entry->t0 + entry->dt > t0))
        return 0;
    return 1;   /* matches */
}

#endif

int XLALCacheSieve(LALCache * cache, INT4 t0, INT4 t1,
                   const char *srcregex, const char *dscregex,
                   const char *urlregex)
{
    UINT4 n = 0;
    UINT4 i;
#ifdef HAVE_REGEX_H
    regex_t srcreg;
    regex_t dscreg;
    regex_t urlreg;
#else /* HAVE_REGEX_H undefined */
    /* can only attempt to match time range */
    if (srcregex || dscregex || urlregex) {
        XLAL_ERROR(XLAL_EFAILED,
                   "Regular expression matching is not supported");
    }
#endif

    if (!cache)
        XLAL_ERROR(XLAL_EFAULT);

#ifdef HAVE_REGEX_H
    if (srcregex)
        regcomp(&srcreg, srcregex, REG_NOSUB);
    if (dscregex)
        regcomp(&dscreg, dscregex, REG_NOSUB);
    if (urlregex)
        regcomp(&urlreg, urlregex, REG_NOSUB);
#endif

    for (i = 0; i < cache->length; ++i) {
        int match;
#ifdef HAVE_REGEX_H
        match =
            XLALCacheEntryMatch(cache->list + i, t0, t1,
                                srcregex ? &srcreg : NULL,
                                dscregex ? &dscreg : NULL,
                                urlregex ? &urlreg : NULL);
#else
        match = XLALCacheEntryMatchTime(cache->list + i, t0, t1);
#endif
        if (match) {
            cache->list[n++] = cache->list[i];
        } else {
            XLALFree(cache->list[i].src);
            XLALFree(cache->list[i].dsc);
            XLALFree(cache->list[i].url);
        }
    }

#ifdef HAVE_REGEX_H
    if (srcregex)
        regfree(&srcreg);
    if (dscregex)
        regfree(&dscreg);
    if (urlregex)
        regfree(&urlreg);
#endif

    if (cache->length != n) {
        cache->list = XLALRealloc(cache->list, n * sizeof(*cache->list));
        if (n && !cache->list)
            XLAL_ERROR(XLAL_EFUNC);
        cache->length = n;
        if (!n)
            XLAL_PRINT_WARNING("No matching entries - zero-length cache");
    }
    return 0;
}


/* 
 * Binary search compare function: to locate the first entry containing the
 * wanted time, or the first file after the requested time if the time is
 * before the beginning or in a gap.
 */
struct LALCacheEntryBsearchKey {
    const void *first;
    double t;
};
static int XLALCacheEntryBsearchCompare(const void *key, const void *ptr)
{
    const struct LALCacheEntryBsearchKey *k = key;
    const LALCacheEntry *entry = ptr;
    double t = k->t;
    /* check if requested time is after this entry */
    if (t > entry->t0 + entry->dt)
        return 1;
    /* check if requested time is before this entry */
    if (t < entry->t0) {
        /* now we need to check:
         * 1. is this the first file?
         * 2. is t after the previous file, i.e., in a gap?
         * if either of these, return 0 */
        if (k->first == ptr || t > entry[-1].t0 + entry[-1].dt)
            return 0;
        return -1;
    }
    /* the requested time is within the bounds of this entry;
     * but is it in the previous file too? */
    if (k->first != ptr && t < entry[-1].t0 + entry[-1].dt)
        return -1;      /* yes it is, so this entry is _after_ the wanted entry */
    return 0;
}

/* 
 * Finds the first entry that contains the requested time, or the first entry
 * after the time if the time is in a gap or before the first entry.  Returns
 * NULL if the time is after the last entry.
 */
LALCacheEntry *XLALCacheEntrySeek(const LALCache * cache, double t)
{
    struct LALCacheEntryBsearchKey key;
    key.first = cache->list;
    key.t = t;
    return bsearch(&key, cache->list, cache->length, sizeof(*cache->list),
                   XLALCacheEntryBsearchCompare);
}

LALFILE *XLALCacheEntryOpen(const LALCacheEntry * entry)
{
    char *nextslash;
    char *nextcolon;
    char *filename;
    LALFILE *fp;
    if (!entry)
        XLAL_ERROR_NULL(XLAL_EFAULT);
    filename = entry->url;
    nextslash = strchr(filename, '/');
    nextcolon = strchr(filename, ':');
    if (nextslash && nextcolon && nextcolon < nextslash
        && 0 == strncmp(nextcolon, "://", 3)) {
        if (strncmp(filename, "file", nextcolon - filename))
            XLAL_ERROR_NULL(XLAL_EIO,
                            "Unsupported protocol in URL %s (only file supported)",
                            entry->url);
        filename = nextcolon + 3;
        nextslash = strchr(filename, '/');
        if (nextslash != filename
            && strncmp(filename, "localhost", nextslash - filename))
            XLAL_ERROR_NULL(XLAL_EIO,
                            "Cannot read files from remote host in URL %s (only localhost supported)",
                            entry->url);
        filename = nextslash;
    }
    fp = XLALFileOpen(filename, "r");
    if (!fp)
        XLAL_ERROR_NULL(XLAL_EIO, "Could not open file %s for input",
                        filename);
    return fp;
}
