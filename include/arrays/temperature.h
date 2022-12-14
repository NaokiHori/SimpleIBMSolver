#if !defined(ARRAYS_TEMPERATURE_H)
#define ARRAYS_TEMPERATURE_H

/* This file is automatically generated by define_array.py */

// check NDIMS (number of spatial dimensions)
#if !defined(NDIMS)
#error "define NDIMS: -DNDIMS=2 or -DNDIMS=3"
#endif

#include "domain.h"

/*** temp ***/
#if NDIMS == 2
#define TEMP_SIZE_0 (isize+2)
#define TEMP_SIZE_1 (jsize+2)
#define TEMP(I, J) (temp[ (J  ) * TEMP_SIZE_0 + (I  ) ])
#endif // NDIMS == 2
#if NDIMS == 3
#define TEMP_SIZE_0 (isize+2)
#define TEMP_SIZE_1 (jsize+2)
#define TEMP_SIZE_2 (ksize+2)
#define TEMP(I, J, K) (temp[ (K  ) * TEMP_SIZE_1 * TEMP_SIZE_0 + (J  ) * TEMP_SIZE_0 + (I  ) ])
#endif // NDIMS == 3
/*** temp ***/

/*** tempforcex ***/
#if NDIMS == 2
#define TEMPFORCEX_SIZE_0 (isize+1)
#define TEMPFORCEX_SIZE_1 (jsize+2)
#define TEMPFORCEX(I, J) (tempforcex[ (J  ) * TEMPFORCEX_SIZE_0 + (I-1) ])
#endif // NDIMS == 2
#if NDIMS == 3
#define TEMPFORCEX_SIZE_0 (isize+1)
#define TEMPFORCEX_SIZE_1 (jsize+2)
#define TEMPFORCEX_SIZE_2 (ksize+2)
#define TEMPFORCEX(I, J, K) (tempforcex[ (K  ) * TEMPFORCEX_SIZE_1 * TEMPFORCEX_SIZE_0 + (J  ) * TEMPFORCEX_SIZE_0 + (I-1) ])
#endif // NDIMS == 3
/*** tempforcex ***/

/*** srctempa ***/
#if NDIMS == 2
#define SRCTEMPA_SIZE_0 (isize+2)
#define SRCTEMPA_SIZE_1 (jsize+2)
#define SRCTEMPA(I, J) (srctempa[ (J  ) * SRCTEMPA_SIZE_0 + (I  ) ])
#endif // NDIMS == 2
#if NDIMS == 3
#define SRCTEMPA_SIZE_0 (isize+2)
#define SRCTEMPA_SIZE_1 (jsize+2)
#define SRCTEMPA_SIZE_2 (ksize+2)
#define SRCTEMPA(I, J, K) (srctempa[ (K  ) * SRCTEMPA_SIZE_1 * SRCTEMPA_SIZE_0 + (J  ) * SRCTEMPA_SIZE_0 + (I  ) ])
#endif // NDIMS == 3
/*** srctempa ***/

/*** srctempb ***/
#if NDIMS == 2
#define SRCTEMPB_SIZE_0 (isize+2)
#define SRCTEMPB_SIZE_1 (jsize+2)
#define SRCTEMPB(I, J) (srctempb[ (J  ) * SRCTEMPB_SIZE_0 + (I  ) ])
#endif // NDIMS == 2
#if NDIMS == 3
#define SRCTEMPB_SIZE_0 (isize+2)
#define SRCTEMPB_SIZE_1 (jsize+2)
#define SRCTEMPB_SIZE_2 (ksize+2)
#define SRCTEMPB(I, J, K) (srctempb[ (K  ) * SRCTEMPB_SIZE_1 * SRCTEMPB_SIZE_0 + (J  ) * SRCTEMPB_SIZE_0 + (I  ) ])
#endif // NDIMS == 3
/*** srctempb ***/

/*** srctempg ***/
#if NDIMS == 2
#define SRCTEMPG_SIZE_0 (isize+2)
#define SRCTEMPG_SIZE_1 (jsize+2)
#define SRCTEMPG(I, J) (srctempg[ (J  ) * SRCTEMPG_SIZE_0 + (I  ) ])
#endif // NDIMS == 2
#if NDIMS == 3
#define SRCTEMPG_SIZE_0 (isize+2)
#define SRCTEMPG_SIZE_1 (jsize+2)
#define SRCTEMPG_SIZE_2 (ksize+2)
#define SRCTEMPG(I, J, K) (srctempg[ (K  ) * SRCTEMPG_SIZE_1 * SRCTEMPG_SIZE_0 + (J  ) * SRCTEMPG_SIZE_0 + (I  ) ])
#endif // NDIMS == 3
/*** srctempg ***/

#endif // ARRAYS_{name}_H
