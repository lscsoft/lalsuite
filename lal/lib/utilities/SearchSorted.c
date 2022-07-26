#include <stddef.h>
#include <sys/types.h>
#include <lal/LALStddef.h>
#include <lal/LALStdlib.h>
#include <lal/Sort.h>

static _LAL_INLINE_ int lt(void *params, const void *a, const void *b, int (*compar)(void *, const void *, const void *))
{
    return compar(params, a, b) < 0;
}

static _LAL_INLINE_ int lte(void *params, const void *a, const void *b, int (*compar)(void *, const void *, const void *))
{
    return compar(params, a, b) <= 0;
}

int XLALIsSorted(void *base, size_t nobj, size_t size, void *params, int (*compar)(void *, const void *, const void *))
{
    char *p = base;
    XLAL_CHECK(base, XLAL_EFAULT);
    XLAL_CHECK(compar, XLAL_EFAULT);
    XLAL_CHECK((ssize_t)size > 0, XLAL_EINVAL);
    while (nobj-- > 1)
        if (compar(params, p + size, p) < 0)
            return 0;
        else
            p += size;
    return 1;
}

ssize_t XLALSearchSorted(const void *key, const void *base, size_t nobj, size_t size, void *params, int (*compar)(void *, const void *, const void *), int side)
{
    size_t imin = 0;
    size_t imax = nobj;
    int (*cmp)(void *, const void *, const void *, int (*)(void *, const void *, const void *)) = side > 0 ? lte : lt;
    XLAL_CHECK(base, XLAL_EFAULT);
    XLAL_CHECK(compar, XLAL_EFAULT);
    XLAL_CHECK((ssize_t)size > 0, XLAL_EINVAL);
    while (imin < imax) {
        size_t i = imin + ((imax - imin) >> 1);
        if (cmp(params, (const char *)base + i * size, key, compar))
            imin = i + 1;
        else
            imax = i;
    }
    if (side || (imin < nobj && compar(params, (const char *)base + imin * size, key) == 0))
        return imin;
    return -1; /* entry not found */
}
