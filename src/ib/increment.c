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
    const double ux1 = p->ux + p->dux;
    const double uy1 = p->uy + p->duy;
    p->dx = 0.5 * dt * (ux0 + ux1);
    p->dy = 0.5 * dt * (uy0 + uy1);
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

