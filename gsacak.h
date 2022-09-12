/*
 * gsacak(s, SA, NULL, NULL, n) //computes only SA
 * gsacak(s, SA, LCP,  NULL, n) //computes SA and LCP
 * gsacak(s, SA, NULL, DA, n)   //computes SA and DA
 * gsacak(s, SA, LCP,  DA, n)   //computes SA, LCP and DA
 *
 */

#ifndef GSACAK_H
#define GSACAK_H

#include <inttypes.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define max(a, b) ((a) > (b) ? (a) : (b))

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef M64
#define M64 0
#endif

typedef int32_t int_t;
typedef uint32_t uint_t;
#define PRIdN PRId32
#define U_MAX UINT32_MAX
#define I_MAX INT32_MAX
#define I_MIN INT32_MIN

/*! @option type of s[0,n-1] for integer alphabets
 *
 *  @constraint sizeof(int_t) >= sizeof(int_text)
 */
typedef uint32_t int_text; // 4N bytes for s[0..n-1]
#define PRIdT PRIu32

#endif
