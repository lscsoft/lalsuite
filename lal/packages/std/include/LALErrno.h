/* <lalVerbatim file="LALErrnoHV">
Author: Cannon, K. C.
Revision: $Id$
</lalVerbatim>
 */

#ifndef _LALERRNO_H
#define _LALERRNO_H

/*
 * Error codes and corresponding error messages.
 */

#define LAL_FAIL_ERR	1
#define LAL_FAIL_MSG	"operation failed"
#define LAL_NULL_ERR	2
#define LAL_NULL_MSG	"unexpected NULL pointer"
#define LAL_NNULL_ERR	3
#define LAL_NNULL_MSG	"unexpected non-NULL pointer"
#define LAL_NOMEM_ERR	4
#define LAL_NOMEM_MSG	"out of memory"
#define LAL_RANGE_ERR	5
#define LAL_RANGE_MSG	"parameter out of range"
#define LAL_BADPARM_ERR 6
#define LAL_BADPARM_MSG "invalid parameter value"

#endif /* _LALERRNO_H */
