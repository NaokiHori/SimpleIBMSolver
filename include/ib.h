#if !defined(IB_H)
#define IB_H

#include "array.h"
#include "domain.h"

typedef struct {
  // radius and density
  double r, d;
  // mass and moment of inertia
  double m, mi;
  // position and displacement
  double x, dx;
  double y, dy;
#if NDIMS == 3
  double z, dz;
#endif
  // translational velocity, surface force, internal inertia, collision force
  double ux, dux, fux, iux[2], cfx[2];
  double uy, duy, fuy, iuy[2], cfy[2];
#if NDIMS == 3
  double uz, duz, fuz, iuz[2], cfz[2];
#endif
  // angular velocity, torque, internal angular inertia
#if NDIMS == 3
  double vx, dvx, tvx, ivx[2];
  double vy, dvy, tvy, ivy[2];
#endif
  double vz, dvz, tvz, ivz[2];
} particle_t;

typedef struct {
#if NDIMS == 2
  // circular objects
#else
  // spherical objects
#endif
  size_t nitems;
  particle_t * particles;
  // response to the momentum field
  array_t ibfx;
  array_t ibfy;
  array_t ibfz;
} ib_t;

#endif // IB_H
