#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "ib.h"
#include "arrays/fluid.h"
#include "arrays/ib.h"



static int kernel_exchange_momentum(const domain_t *domain, const fluid_t *fluid, const double dt, ib_t *ib){
  // \int -fx dV
  // \int -fy dV
  // \int ( + ry fx - rx fy ) dV
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ioffs = domain->offsets[0];
  const int joffs = domain->offsets[1];
  const double ly = domain->lengths[1];
  const double dx = domain->dx;
  const double dy = domain->dy;
  const double xoffs = dx * ioffs;
  const double yoffs = dy * joffs;
  const double *ux = fluid->ux;
  const double *uy = fluid->uy;
  const int n_particles = ib->n_particles;
  double *dux = ib->dux;
  double *duy = ib->duy;
  memset(dux, 0, DUX_SIZE_0 * DUX_SIZE_1 * sizeof(double));
  memset(duy, 0, DUY_SIZE_0 * DUY_SIZE_1 * sizeof(double));
  for(int n = 0; n < n_particles; n++){
    particle_t *p = ib->particles[n];
    // constant parameters
    const double pr  = p->r;
    const double px  = p->x;
    const double py  = p->y;
    const double pux = p->ux;
    const double puy = p->uy;
    const double pvz = p->vz;
    const double pm  = ib_compute_mass(1., pr);
    const double pim = ib_compute_moment_of_inertia(1., pr);
    // buffers (for simplicity)
    double fux = 0.;
    double fuy = 0.;
    double tvz = 0.;
    for(int periodic = -1; periodic <= 1; periodic++){
      double py_ = py + ly * periodic;
      int imin, imax, jmin, jmax;
      ib_decide_loop_size(1, isize, dx, pr, px -xoffs, &imin, &imax);
      ib_decide_loop_size(1, jsize, dy, pr, py_-yoffs, &jmin, &jmax);
      for(int j = jmin; j <= jmax; j++){
        double y = yoffs + 0.5 * (2 * j - 1) * dy;
        for(int i = imin; i <= imax; i++){
          double x = xoffs + 0.5 * (2 * i - 1) * dx;
          double w = ib_s_weight(fmin(dx, dy), pr, px, py_, x, y);
          double ux_p = pux-pvz*(y-py_);
          double uy_p = puy+pvz*(x-px );
          double ux_f = 0.5*(UX(i  , j  )+UX(i+1, j  ));
          double uy_f = 0.5*(UY(i  , j  )+UY(i  , j+1));
          double fx = w*(ux_p-ux_f)/dt;
          double fy = w*(uy_p-uy_f)/dt;
          DUX(i, j) += fx*dt;
          DUY(i, j) += fy*dt;
          fux -= +(fx*dx*dy)/pm;
          tvz -= -(y-py_)*(fx*dx*dy)/pim;
          fuy -= +(fy*dx*dy)/pm;
          tvz -= +(x-px)*(fy*dx*dy)/pim;
        }
      }
    }
    fluid_update_boundaries_p(domain, dux);
    fluid_update_boundaries_p(domain, duy);
    // assign results
    p->fux = fux;
    p->fuy = fuy;
    p->tvz = tvz;
  }
  return 0;
}

static int synchronise_information(ib_t *ib){
  const int n_particles = ib->n_particles;
  particle_t **particles = ib->particles;
  // prepare message buffer
  double *buf = common_calloc(3*n_particles, sizeof(double));
  // pack
  for(int n = 0; n < n_particles; n++){
    particle_t *p = particles[n];
    buf[3*n+0] = p->fux;
    buf[3*n+1] = p->fuy;
    buf[3*n+2] = p->tvz;
  }
  // sum up all
  MPI_Allreduce(MPI_IN_PLACE, buf, 3*n_particles, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // unpack
  for(int n = 0; n < n_particles; n++){
    particle_t *p = particles[n];
    p->fux = buf[3*n+0];
    p->fuy = buf[3*n+1];
    p->tvz = buf[3*n+2];
  }
  common_free(buf);
  return 0;
}


int ib_exchange_momentum(const domain_t *domain, const fluid_t *fluid, const double dt, ib_t *ib){
  // exchange (translational and angular) momenta with each particle
  kernel_exchange_momentum(domain, fluid, dt, ib);
  // communicate updated information
  synchronise_information(ib);
  return 0;
}

