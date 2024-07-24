#include "memory.h"
#include "domain.h"
#include "fluid.h"
#include "ib.h"
#include "ib_solver.h"
#include "internal.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"

static int compute_inertia(
    const domain_t * domain,
    const size_t index,
    const fluid_t * fluid,
    ib_t * ib
){
  // \int ux dV
  // \int uy dV
  // \int ( - ry ux + rx uy ) dV
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
  for(size_t n = 0; n < nitems; n++){
    particle_t * p = particles + n;
    // constant parameters
    const double pr  = p->r;
    const double px  = p->x + p->dx;
    const double py  = p->y + p->dy;
    const double pm  = p->m;
    const double pmi = p->mi;
    // targets
    double * iux = &p->iux[index];
    double * iuy = &p->iuy[index];
    double * ivz = &p->ivz[index];
    *iux = 0.;
    *iuy = 0.;
    *ivz = 0.;
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
          const double w = ib_v_weight(
              deltas,
              pr,
              (double [NDIMS]){px, py_},
              (double [NDIMS]){x, y}
          );
          const double velx = + 0.5 * UX(i  , j  ) + 0.5 * UX(i+1, j  );
          const double vely = + 0.5 * UY(i  , j  ) + 0.5 * UY(i  , j+1);
          const double valx = w * velx * cellsize;
          const double valy = w * vely * cellsize;
          *iux += valx / pm;
          *iuy += valy / pm;
          *ivz += - (y - py_) * valx / pmi;
          *ivz += + (x - px ) * valy / pmi;
        }
      }
    }
  }
  return 0;
}

static int synchronise_information(
    const domain_t * domain,
    const size_t index,
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
    buf[nvars * n + 0] = p->iux[index];
    buf[nvars * n + 1] = p->iuy[index];
    buf[nvars * n + 2] = p->ivz[index];
  }
  // sum up all
  MPI_Allreduce(MPI_IN_PLACE, buf, nvars * nitems, MPI_DOUBLE, MPI_SUM, comm_cart);
  // unpack
  for(size_t n = 0; n < nitems; n++){
    particle_t * p = particles + n;
    p->iux[index] = buf[nvars * n + 0];
    p->iuy[index] = buf[nvars * n + 1];
    p->ivz[index] = buf[nvars * n + 2];
  }
  memory_free(buf);
  return 0;
}

int ib_compute_inertia(
    const domain_t * domain,
    const size_t index,
    const fluid_t * fluid,
    ib_t * ib
){
  if(0 != compute_inertia(domain, index, fluid, ib)){
    return 1;
  }
  if(0 != synchronise_information(domain, index, ib)){
    return 1;
  }
  return 0;
}

