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

/******************************************************************************/

/** @brief computes the suffix array of string s[0..n-1]
 *
 *  @param k	alphabet size
 */
int sacak_int(int_text *s, uint_t *SA, uint_t n, uint_t k);

/******************************************************************************/

/** @brief Computes the suffix array SA (LCP, DA) of T^cat in s[0..n-1]
 *
 *  @param s		input concatenated string, using separators s[i]=1 and
 * with s[n-1]=0
 *  @param SA		suffix array
 *  @param LCP	LCP array
 *  @param DA		Document array
 *  @param n	string length
 *
 *  @return depth of the recursive calls.
 */
int gsacak(unsigned char *s, uint_t *SA, int_t *LCP, int_t *DA, uint_t n);

/** @brief Computes the suffix array SA (LCP, DA) of T^cat in s[0..n-1]
 *
 *  @param s    input concatenated string, with s[n-1]=0
 *  @param K	alphabet size
 */
int gsacak_int(
    int_text *s, uint_t *SA, int_t *LCP, int_t *DA, uint_t n, uint_t k
);

/******************************************************************************/

int_t gSACA_K(
    uint_t *s, uint_t *SA, uint_t n, unsigned int K, int cs, uint_t separator,
    int level
);

int_t SACA_K(
    int_t *s, uint_t *SA, uint_t n, unsigned int K, uint_t m, int cs, int level
);

#endif
