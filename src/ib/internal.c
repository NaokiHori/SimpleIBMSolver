#include <math.h>
#include "ib.h"
#include "internal.h"

// sharpness of diffusive surface
static const double beta = 2.;

// number of grid points outside circles
static const int nextra = 3;

static double min(
    const double vals[NDIMS]
){
  double val = 1.e+16;
  for(size_t dim = 0; dim < NDIMS; dim++){
    val = fmin(val, vals[dim]);
  }
  return val;
}

static double compute_signed_dist(
    const double delta,
    const double radius,
    const double center[NDIMS],
    const double point[NDIMS]
){
  const double deltas[NDIMS] = {
    point[0] - center[0],
    point[1] - center[1],
#if NDIMS == 3
    point[2] - center[2],
#endif
  };
  double dist = sqrt(
      + deltas[0] * deltas[0]
      + deltas[1] * deltas[1]
#if NDIMS == 3
      + deltas[2] * deltas[2]
#endif
  );
  // convert it such that inside: >0, outside: <0
  dist = radius - dist;
  return dist / delta;
}

double ib_s_weight(
    const double deltas[NDIMS],
    const double pr,
    const double center[NDIMS],
    const double point[NDIMS]
){
  const double dist = compute_signed_dist(min(deltas), pr, center, point);
  const double tanh_ = tanh(beta * dist);
  return 0.5 * beta * (1. - tanh_ * tanh_);
}

double ib_v_weight(
    const double deltas[NDIMS],
    const double pr,
    const double center[NDIMS],
    const double point[NDIMS]
){
  const double dist = compute_signed_dist(min(deltas), pr, center, point);
  const double tanh_ = tanh(beta * dist);
  return 0.5 * (1. + tanh_);
}

int ib_decide_loop_size(
    const int lbound,
    const int ubound,
    const double delta,
    const double radius,
    const double center,
    int * min,
    int * max
){
  *min = (int)( floor((center - radius) / delta) - nextra );
  *max = (int)( ceil ((center + radius) / delta) + nextra );
  *min = *min < lbound ? lbound : *min;
  *max = *max > ubound ? ubound : *max;
  return 0;
}

