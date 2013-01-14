/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Robert Adam Mercer
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

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALFrameU.h>

#ifndef PATH_MAX
#define PATH_MAX FILENAME_MAX
#endif

#define H0 2.3e-18      /* Hubble constant, s */
#define SITE_MAX 26     /* in our dreams... */

#define FAILURE(...) do { fprintf(stderr, __VA_ARGS__); exit(1); } while (0)

int charcmp(const void *c1, const void *c2);
char *fnametodsc(const char *fname);

int main(int argc, char *argv[])
{
    char stdio[] = "-";
    char *defaultfilev[1] = { stdio };
    int filec = (argc == 1 ? 1 : argc - 1);
    char **filev = (argc == 1 ? defaultfilev : &argv[1]);
    int f;

    for (f = 0; f < filec; ++f) {
        char *fname = filev[f];
        char *dsc;
        LALFrameUFrFile *frfile;
        LALFrameUFrTOC *toc;
        char sites[SITE_MAX + 1] = "";
        char path[PATH_MAX + 1] = "";
        size_t ndet;
        size_t nadc;
        size_t nsim;
        size_t nproc;
        size_t nframe;
        size_t pos;
        size_t det;
        double t0min = +1.0 / H0;       /* ridiculously far in future */
        double t1max = -1.0 / H0;       /* ridiculously far in past */
        int dt;

        frfile = XLALFrameUFrFileOpen(fname, "r");
        if (!frfile)
            FAILURE("file %s not found\n", fname);

        toc = XLALFrameUFrTOCRead(frfile);
        if (!toc)
            FAILURE("no TOC found\n");


        nframe = XLALFrameUFrTOCQueryNFrame(toc);
        if ((int) (nframe) <= 0)
            FAILURE("no frames found\n");

        ndet = XLALFrameUFrTOCQueryDetectorN(toc);
        nadc = XLALFrameUFrTOCQueryAdcN(toc);
        nsim = XLALFrameUFrTOCQuerySimN(toc);
        nproc = XLALFrameUFrTOCQueryProcN(toc);

        for (det = 0; det < ndet; ++det) {
            LALFrameUFrDetector *detector;
            const char *prefix;
            const char *name;
            name = XLALFrameUFrTOCQueryDetectorName(toc, det);
            detector = XLALFrameUFrDetectorRead(frfile, name);
            prefix = XLALFrameUFrDetectorQueryPrefix(detector);
            /* add site if it is new */
            if (prefix && isupper(*prefix))
                if (strchr(sites, *prefix) == NULL)
                    strncat(sites, prefix, 1);
            XLALFrameUFrDetectorFree(detector);
        }
        /* sort letters in site string */
        if (*sites)
            qsort(sites, strlen(sites), sizeof(*sites), charcmp);
        else
            strcpy(sites, "-");

        for (pos = 0; pos < nframe; ++pos) {
            LALFrameUFrChan *channel;
            double tip;
            double tfp;
            double t0;
            double t1;
            size_t chan;
            tfp = XLALFrameUFrTOCQueryGTimeModf(&tip, toc, pos);
            t0 = tip + tfp;
            t1 = t0 + XLALFrameUFrTOCQueryDt(toc, pos);
            if (t0 < t0min)
                t0min = t0;
            if (t1 > t1max)
                t1max = t1;

            for (chan = 0; chan < nadc; ++chan) {
                const char *name;
                double toff;
                name = XLALFrameUFrTOCQueryAdcName(toc, chan);
                channel = XLALFrameUFrChanRead(frfile, name, pos);
                if (!channel)
                    continue;
                toff = XLALFrameUFrChanQueryTimeOffset(channel);
                if (t0 + toff < t0min)
                    t0min = t0 + toff;
                if (t1 + toff > t1max)
                    t1max = t1 + toff;
                XLALFrameUFrChanFree(channel);
            }

            for (chan = 0; chan < nsim; ++chan) {
                const char *name;
                double toff;
                name = XLALFrameUFrTOCQuerySimName(toc, chan);
                channel = XLALFrameUFrChanRead(frfile, name, pos);
                if (!channel)
                    continue;
                toff = XLALFrameUFrChanQueryTimeOffset(channel);
                if (t0 + toff < t0min)
                    t0min = t0 + toff;
                if (t1 + toff > t1max)
                    t1max = t1 + toff;
                XLALFrameUFrChanFree(channel);
            }

            for (chan = 0; chan < nproc; ++chan) {
                const char *name;
                double toff;
                name = XLALFrameUFrTOCQueryProcName(toc, chan);
                channel = XLALFrameUFrChanRead(frfile, name, pos);
                if (!channel)
                    continue;
                toff = XLALFrameUFrChanQueryTimeOffset(channel);
                if (t0 + toff < t0min)
                    t0min = t0 + toff;
                if (t1 + toff > t1max)
                    t1max = t1 + toff;
                XLALFrameUFrChanFree(channel);
            }
        }

        dsc = fnametodsc(fname);
        if (!realpath(fname, path))
            strcpy(path, "-");
        dt = ceil(t1max) - floor(t0min);
        printf("%s", sites);
        printf("\t%s", dsc);
        printf("\t%d", (int) floor(t0min));
        printf("\t%d", dt > 0 ? dt : 1);
        if (strcmp(path, "-") == 0)
            printf("\t-");
        else
            printf("\tfile://%s", path);
        printf("\n");

        free(dsc);
        XLALFrameUFrTOCFree(toc);
        XLALFrameUFrFileClose(frfile);
    }

    return 0;
}

int charcmp(const void *c1, const void *c2)
{
    char a = *((const char *) c1);
    char b = *((const char *) c2);
    return (a > b) - (a < b);
}

char *fnametodsc(const char *fname)
{
    if (fname && *fname && strcmp(fname, "-") && strlen(fname) < FILENAME_MAX) {
        char src[FILENAME_MAX];
        char dsc[FILENAME_MAX];
        const char *base;
        char *s;
        unsigned int t0;
        unsigned int dt;
        int n;

        /* get basename */
        base = strrchr(fname, '/');
        base = base ? base + 1 : fname;

        /* see if basename is in canonical form */
        n = sscanf(base, "%[A-Z]-%[a-zA-Z0-9_+#]-%u-%u.gwf", src, dsc, &t0, &dt);
        if (n == 4)     /* correct number of conversions */
            return strdup(dsc); /* return description field */

        /* try to convert the full basename */
        base = strrchr(fname, '/');
        base = base ? base + 1 : fname;

        /* translate spaces & hyphens to octothorpes & underscores */
        s = dsc;
        while (*base && *base != '.')
            switch (*base) {
            case ' ':
                *s++ = '#';
                ++base;
                break;
            case '-':
                *s++ = '_';
                ++base;
                break;
            default:
                *s++ = *base++;
                break;
            }
        *s = '\0';
        if (strlen(dsc))
            return strdup(dsc);
    }

    /* all else fails... */
    return strdup("-");
}
