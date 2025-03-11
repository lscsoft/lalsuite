#include <stddef.h>
#include <sys/types.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/Sort.h>

static void msort(void *base, size_t nobj, size_t size, void *params, int (*compar)(void *, const void *, const void *), void *work)
{
    if (nobj > 1) {
        typedef char type[size];
        size_t n1 = nobj >> 1;
        size_t n2 = nobj - n1;
        type *a1 = (type *)base;
        type *a2 = (type *)base + n1;
        type *tmp = (type *)work;

        msort(a1, n1, size, params, compar, work);
        msort(a2, n2, size, params, compar, work);

        while (n1 > 0 && n2 > 0) {
            if (compar(params, a1, a2) > 0) {
                memcpy(tmp++, a2++, size);
                --n2;
            } else {
                memcpy(tmp++, a1++, size);
                --n1;
            }
        }
        if (n1 > 0)
            memcpy(tmp, a1, n1 * size);
        memcpy(base, work, (nobj - n2) * size);
    }
    return;
}

int XLALMergeSort(void *base, size_t nobj, size_t size, void *params, int (*compar)(void *, const void *, const void *))
{
    void *work;
    if (nobj < 2) /* 0 or 1 objects are already sorted */
        return 0;
    XLAL_CHECK(base, XLAL_EFAULT);
    XLAL_CHECK(compar, XLAL_EFAULT);
    XLAL_CHECK((ssize_t)size > 0, XLAL_EINVAL);
    work = LALMalloc(nobj * size);
    XLAL_CHECK(work, XLAL_ENOMEM);
    msort(base, nobj, size, params, compar, work);
    LALFree(work);
    return 0;
}
