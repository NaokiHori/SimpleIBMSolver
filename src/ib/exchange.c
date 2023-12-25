#include <math.h>
#include <mpi.h>
#include "runge_kutta.h"
#include "memory.h"
#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "ib.h"
#include "ib_solver.h"
#include "internal.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/ib/ibfx.h"
#include "array_macros/ib/ibfy.h"

static int update_boundaries(
    const domain_t * domain,
    array_t * arr
){
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, arr)){
    return 1;
  }
  return 0;
}

static int exchange_momentum(
    const domain_t * domain,
    const fluid_t * fluid,
    const double dt,
    ib_t * ib
){
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
  const double deltas[NDIMS] = {dx, dy};
  const double cellsize = dx * dy;
  const double xoffs = dx * ioffs;
  const double yoffs = dy * joffs;
  const double * ux = fluid->ux.data;
  const double * uy = fluid->uy.data;
  const size_t nitems = ib->nitems;
  particle_t * particles = ib->particles;
  double * ibfx = ib->ibfx.data;
  double * ibfy = ib->ibfy.data;
  for(size_t n = 0; n < nitems; n++){
    particle_t * p = particles + n;
    // constant parameters
    const double pr  = p->r;
    const double px  = p->x;
    const double py  = p->y;
    const double pux = p->ux;
    const double puy = p->uy;
    const double pvz = p->vz;
    // mass and moment of inertia
    const double pm  = p->m;
    const double pmi = p->mi;
    // targets
    double * fux = &p->fux;
    double * fuy = &p->fuy;
    double * tvz = &p->tvz;
    *fux = 0.;
    *fuy = 0.;
    *tvz = 0.;
    for(int yper = -1; yper <= 1; yper++){
      const double py_ = py + ly * yper;
      int imin = 0;
      int imax = 0;
      int jmin = 0;
      int jmax = 0;
      ib_decide_loop_size(1, isize, dx, pr, px  - xoffs, &imin, &imax);
      ib_decide_loop_size(1, jsize, dy, pr, py_ - yoffs, &jmin, &jmax);
      for(int j = jmin; j <= jmax; j++){
        const double y = yoffs + 0.5 * (2 * j - 1) * dy;
        for(int i = imin; i <= imax; i++){
          const double x = xoffs + 0.5 * (2 * i - 1) * dx;
          const double w = ib_s_weight(
              deltas,
              pr,
              (double [NDIMS]){px, py_},
              (double [NDIMS]){x, y}
          );
          const double ux_p = pux - pvz * (y - py_);
          const double uy_p = puy + pvz * (x - px );
          const double ux_f = + 0.5 * UX(i  , j  ) + 0.5 * UX(i+1, j  );
          const double uy_f = + 0.5 * UY(i  , j  ) + 0.5 * UY(i  , j+1);
          const double fx = w * (ux_p - ux_f) / dt;
          const double fy = w * (uy_p - uy_f) / dt;
          IBFX(i, j) += fx * dt;
          IBFY(i, j) += fy * dt;
          *fux -= + fx * cellsize / pm;
          *fuy -= + fy * cellsize / pm;
          *tvz -= - (y - py_) * (fx * cellsize) / pmi;
          *tvz -= + (x - px ) * (fy * cellsize) / pmi;
        }
      }
    }
    update_boundaries(domain, &ib->ibfx);
    update_boundaries(domain, &ib->ibfy);
  }
  return 0;
}

static int update_predicted_velocity(
    const domain_t * domain,
    fluid_t * fluid,
    const ib_t * ib
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  double * ux = fluid->ux.data;
  double * uy = fluid->uy.data;
  const double * ibfx = ib->ibfx.data;
  const double * ibfy = ib->ibfy.data;
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      UX(i, j) +=
        + 0.5 * IBFX(i-1, j  )
        + 0.5 * IBFX(i  , j  );
    }
  }
  if(0 != fluid_update_boundaries_ux(domain, &fluid->ux)){
    return 1;
  }
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      UY(i, j) +=
        + 0.5 * IBFY(i  , j-1)
        + 0.5 * IBFY(i  , j  );
    }
  }
  if(0 != fluid_update_boundaries_uy(domain, &fluid->uy)){
    return 1;
  }
  return 0;
}

static int synchronise_information(
    const domain_t * domain,
    ib_t * ib
){
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  const size_t nitems = ib->nitems;
  particle_t * particles = ib->particles;
  // prepare message buffer
  const size_t nvars = 3;
  double * buf = memory_calloc(nvars * nitems, sizeof(double));
  // pack
  for(size_t n = 0; n < nitems; n++){
    particle_t * p = particles + n;
    buf[nvars * n + 0] = p->fux;
    buf[nvars * n + 1] = p->fuy;
    buf[nvars * n + 2] = p->tvz;
  }
  // sum up all
  MPI_Allreduce(MPI_IN_PLACE, buf, nvars * nitems, MPI_DOUBLE, MPI_SUM, comm_cart);
  // unpack
  for(size_t n = 0; n < nitems; n++){
    particle_t * p = particles + n;
    p->fux = buf[nvars * n + 0];
    p->fuy = buf[nvars * n + 1];
    p->tvz = buf[nvars * n + 2];
  }
  memory_free(buf);
  return 0;
}

int ib_exchange_momentum(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid,
    ib_t * ib
){
  if(0 != exchange_momentum(domain, fluid, rkcoefs[rkstep][rk_g] * dt, ib)){
    return 1;
  }
  if(0 != update_predicted_velocity(domain, fluid, ib)){
    return 1;
  }
  if(0 != synchronise_information(domain, ib)){
    return 1;
  }
  return 0;
}

