#if !defined(IB_INTERNAL_H)
#define IB_INTERNAL_H

#include "domain.h"

extern int ib_decide_loop_size(
    const int lbound,
    const int ubound,
    const double delta,
    const double radius,
    const double center,
    int * min,
    int * max
);

extern double ib_s_weight(
    const double deltas[NDIMS],
    const double pr,
    const double center[NDIMS],
    const double point[NDIMS]
);

extern double ib_v_weight(
    const double deltas[NDIMS],
    const double pr,
    const double center[NDIMS],
    const double point[NDIMS]
);

#endif // IB_INTERNAL_H
