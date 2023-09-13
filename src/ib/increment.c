#include "runge_kutta.h"
#include "ib.h"
#include "ib_solver.h"

static int increment_particle_velocities(
    const double dt,
    ib_t * ib
){
  const size_t nitems = ib->nitems;
  particle_t * particles = ib->particles;
  for(size_t n = 0; n < nitems; n++){
    particle_t * p = particles + n;
    // compute new increments
    p->dux =
      + 0.5 * dt * (p->cfx[0] + p->cfx[1]) // collision
      +       dt * p->fux                  // boundary force
      + p->iux[1] - p->iux[0]              // internal inertia
      ;
    p->duy =
      + 0.5 * dt * (p->cfy[0] + p->cfy[1]) // collision
      +       dt * p->fuy                  // boundary force
      + p->iuy[1] - p->iuy[0]              // internal inertia
      ;
#if NDIMS == 3
    p->duz =
      + 0.5 * dt * (p->cfz[0] + p->cfz[1]) // collision
      +       dt * p->fuz                  // boundary force
      + p->iuz[1] - p->iuz[0]              // internal inertia
      ;
#endif
#if NDIMS == 3
    p->dvx =
      +       dt * p->tvx                  // boundary torque
      + p->ivx[1] - p->ivx[0]              // internal inertia
      ;
    p->dvy =
      +       dt * p->tvy                  // boundary torque
      + p->ivy[1] - p->ivy[0]              // internal inertia
      ;
#endif
    p->dvz =
      +       dt * p->tvz                  // boundary torque
      + p->ivz[1] - p->ivz[0]              // internal inertia
      ;
  }
  return 0;
}

static int increment_particle_positions(
    const double dt,
    ib_t * ib
){
  const size_t nitems = ib->nitems;
  particle_t * particles = ib->particles;
  for(size_t n = 0; n < nitems; n++){
    particle_t * p = particles + n;
    const double ux0 = p->ux;
    const double uy0 = p->uy;
#if NDIMS == 3
    const double uz0 = p->uz;
#endif
    const double ux1 = p->ux + p->dux;
    const double uy1 = p->uy + p->duy;
#if NDIMS == 3
    const double uz1 = p->uz + p->duz;
#endif
    p->dx = 0.5 * dt * (ux0 + ux1);
    p->dy = 0.5 * dt * (uy0 + uy1);
#if NDIMS == 3
    p->dz = 0.5 * dt * (uz0 + uz1);
#endif
  }
  return 0;
}

int ib_increment_particles(
    const size_t rkstep,
    double dt,
    ib_t * ib
){
  dt *= rkcoefs[rkstep][rk_g];
  increment_particle_velocities(dt, ib);
  increment_particle_positions(dt, ib);
  return 0;
}

