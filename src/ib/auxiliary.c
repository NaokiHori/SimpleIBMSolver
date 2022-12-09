#include <math.h>
#include "common.h" // M_PI
#include "ib.h"


// number of grid points outside circles
#define NEXTRA 3
// sharpness of diffusive surface
#define BETA 2.

int ib_decide_loop_size(const int lbound, const int ubound, const double grid_size, const double radius, const double grav_center, int *min, int *max){
  *min = (int)(floor((grav_center-radius)/grid_size)-NEXTRA);
  *max = (int)(ceil ((grav_center+radius)/grid_size)+NEXTRA);
  *min = *min < lbound ? lbound : *min;
  *max = *max > ubound ? ubound : *max;
  return 0;
}


static double compute_signed_dist(const double grid_size, const double radius, const double px, const double py, const double pz, const double x, const double y, const double z){
  const double dx = x - px;
  const double dy = y - py;
  const double dz = z - pz;
  // inside: >0, outside: <0
  const double dist
    = radius
    - sqrt(
      + pow(dx, 2.)
      + pow(dy, 2.)
      + pow(dz, 2.)
    );
  return dist/grid_size;
}

double ib_s_weight(const double grid_size, const double radius, const double px, const double py, const double pz, const double x, const double y, const double z){
  double dist = compute_signed_dist(grid_size, radius, px, py, pz, x, y, z);
  return 0.5 * BETA * ( 1. - pow( tanh( BETA * dist ), 2.) );
}

double ib_v_weight(const double grid_size, const double radius, const double px, const double py, const double pz, const double x, const double y, const double z){
  double dist = compute_signed_dist(grid_size, radius, px, py, pz, x, y, z);
  return 0.5 * ( 1. + tanh( BETA * dist ) );
}

double ib_compute_volume(const double radius){
  return 4. / 3. * M_PI * pow(radius, 2.);
}

double ib_compute_mass(const double den, const double radius){
  double vol = ib_compute_volume(radius);
  return den * vol;
}

double ib_compute_moment_of_inertia(const double den, const double radius){
  double mass = ib_compute_mass(den, radius);
  return 0.4 * mass * pow(radius, 2.);
}


#undef NEXTRA
#undef BETA

