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
  // translational velocity, surface force, internal inertia, collision force
  double ux, dux, fux, iux[2], cfx[2];
  double uy, duy, fuy, iuy[2], cfy[2];
  // angular velocity, torque, internal angular inertia
  double vz, dvz, tvz, ivz[2];
} particle_t;

typedef struct {
  // circular objects
  size_t nitems;
  particle_t * particles;
  // response to the momentum field
  array_t ibfx;
  array_t ibfy;
  array_t ibfz;
} ib_t;

#endif // IB_H
