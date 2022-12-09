#include <math.h>
#include <float.h>
#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "ib.h"



static int increment_particle_velocities(const int rkstep, const double dt, ib_t *ib, double *residual){
  const double gamma = RKCOEFS[rkstep].gamma;
  const int n_particles = ib->n_particles;
  for(int n = 0; n < n_particles; n++){
    particle_t *p = ib->particles[n];
    // store previous increments
    double dux_prev = p->dux;
    double duy_prev = p->duy;
    double duz_prev = p->duz;
    double dvx_prev = p->dvx;
    double dvy_prev = p->dvy;
    double dvz_prev = p->dvz;
    // compute new increments
    p->dux =
      +0.5*gamma*dt*(p->cfx[0]+p->cfx[1]) // collision
      +          dt* p->fux               // boundary force
      +(p->iux[1]-p->iux[0])              // internal inertia
      ;
    p->duy =
      +0.5*gamma*dt*(p->cfy[0]+p->cfy[1]) // collision
      +          dt* p->fuy               // boundary force
      +(p->iuy[1]-p->iuy[0])              // internal inertia
      ;
    p->duz =
      +0.5*gamma*dt*(p->cfz[0]+p->cfz[1]) // collision
      +          dt* p->fuz               // boundary force
      +(p->iuz[1]-p->iuz[0])              // internal inertia
      ;
    p->dvx =
      +          dt* p->tvx               // boundary torque
      +(p->ivx[1]-p->ivx[0])              // internal inertia
      ;
    p->dvy =
      +          dt* p->tvy               // boundary torque
      +(p->ivy[1]-p->ivy[0])              // internal inertia
      ;
    p->dvz =
      +          dt* p->tvz               // boundary torque
      +(p->ivz[1]-p->ivz[0])              // internal inertia
      ;
    // check convergence
    *residual = fmax(*residual, fabs(p->dux - dux_prev));
    *residual = fmax(*residual, fabs(p->duy - duy_prev));
    *residual = fmax(*residual, fabs(p->duz - duz_prev));
    *residual = fmax(*residual, fabs(p->dvx - dvx_prev));
    *residual = fmax(*residual, fabs(p->dvy - dvy_prev));
    *residual = fmax(*residual, fabs(p->dvz - dvz_prev));
  }
  return 0;
}

static int increment_particle_positions(const int rkstep, const double dt, ib_t *ib){
  const double gamma = RKCOEFS[rkstep].gamma;
  const int n_particles = ib->n_particles;
  for(int n = 0; n < n_particles; n++){
    particle_t *p = ib->particles[n];
    double ux0 = p->ux;
    double uy0 = p->uy;
    double uz0 = p->uz;
    double ux1 = p->ux+p->dux;
    double uy1 = p->uy+p->duy;
    double uz1 = p->uz+p->duz;
    p->dx = 0.5*gamma*dt*(ux0+ux1);
    p->dy = 0.5*gamma*dt*(uy0+uy1);
    p->dz = 0.5*gamma*dt*(uz0+uz1);
  }
  return 0;
}


int ib_increment_particles(const int rkstep, const double dt, ib_t *ib, double *residual){
  *residual = 0.;
  increment_particle_velocities(rkstep, dt, ib, residual);
  increment_particle_positions(rkstep, dt, ib);
  MPI_Allreduce(MPI_IN_PLACE, residual, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return 0;
}

