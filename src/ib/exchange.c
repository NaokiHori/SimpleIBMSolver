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
#if NDIMS == 3
#include "array_macros/fluid/uz.h"
#endif
#include "array_macros/ib/ibfx.h"
#include "array_macros/ib/ibfy.h"
#if NDIMS == 3
#include "array_macros/ib/ibfz.h"
#endif

static int update_boundaries(
    const domain_t * domain,
    array_t * arr
){
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
#if NDIMS == 3
    MPI_DOUBLE,
#endif
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, arr)){
    return 1;
  }
#if NDIMS == 3
  if(0 != halo_communicate_in_z(domain, dtypes + 1, arr)){
    return 1;
  }
#endif
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
#if NDIMS == 3
  // \int -fz dV
#endif
#if NDIMS == 3
  // \int ( + rz fy - ry fz ) dV
  // \int ( + rx fz - rz fx ) dV
#endif
  // \int ( + ry fx - rx fy ) dV
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const int ioffs = domain->offsets[0];
  const int joffs = domain->offsets[1];
#if NDIMS == 3
  const int koffs = domain->offsets[2];
#endif
  const double ly = domain->lengths[1];
#if NDIMS == 3
  const double lz = domain->lengths[2];
#endif
  const double dx = domain->dx;
  const double dy = domain->dy;
#if NDIMS == 3
  const double dz = domain->dz;
#endif
#if NDIMS == 2
  const double deltas[NDIMS] = {dx, dy};
#else
  const double deltas[NDIMS] = {dx, dy, dz};
#endif
#if NDIMS == 2
  const double cellsize = dx * dy;
#else
  const double cellsize = dx * dy * dz;
#endif
  const double xoffs = dx * ioffs;
  const double yoffs = dy * joffs;
#if NDIMS == 3
  const double zoffs = dz * koffs;
#endif
  const double * ux = fluid->ux.data;
  const double * uy = fluid->uy.data;
#if NDIMS == 3
  const double * uz = fluid->uz.data;
#endif
  const size_t nitems = ib->nitems;
  particle_t * particles = ib->particles;
  double * ibfx = ib->ibfx.data;
  double * ibfy = ib->ibfy.data;
#if NDIMS == 3
  double * ibfz = ib->ibfz.data;
#endif
  for(size_t n = 0; n < nitems; n++){
    particle_t * p = particles + n;
    // constant parameters
    const double pr  = p->r;
    const double px  = p->x;
    const double py  = p->y;
#if NDIMS == 3
    const double pz  = p->z;
#endif
    const double pux = p->ux;
    const double puy = p->uy;
#if NDIMS == 3
    const double puz = p->uz;
#endif
#if NDIMS == 3
    const double pvx = p->vx;
    const double pvy = p->vy;
#endif
    const double pvz = p->vz;
    // mass and moment of inertia
    const double pm  = p->m;
    const double pmi = p->mi;
    // targets
    double * fux = &p->fux;
    double * fuy = &p->fuy;
#if NDIMS == 3
    double * fuz = &p->fuz;
#endif
#if NDIMS == 3
    double * tvx = &p->tvx;
    double * tvy = &p->tvy;
#endif
    double * tvz = &p->tvz;
    *fux = 0.;
    *fuy = 0.;
#if NDIMS == 3
    *fuz = 0.;
#endif
#if NDIMS == 3
    *tvx = 0.;
    *tvy = 0.;
#endif
    *tvz = 0.;
#if NDIMS == 2
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
#else
    for(int zper = -1; zper <= 1; zper++){
      double pz_ = pz + lz * zper;
      for(int yper = -1; yper <= 1; yper++){
        double py_ = py + ly * yper;
        int imin = 0;
        int imax = 0;
        int jmin = 0;
        int jmax = 0;
        int kmin = 0;
        int kmax = 0;
        ib_decide_loop_size(1, isize, dx, pr, px  - xoffs, &imin, &imax);
        ib_decide_loop_size(1, jsize, dy, pr, py_ - yoffs, &jmin, &jmax);
        ib_decide_loop_size(1, ksize, dz, pr, pz_ - zoffs, &kmin, &kmax);
        for(int k = kmin; k <= kmax; k++){
          const double z = zoffs + 0.5 * (2 * k - 1) * dz;
          for(int j = jmin; j <= jmax; j++){
            const double y = yoffs + 0.5 * (2 * j - 1) * dy;
            for(int i = imin; i <= imax; i++){
              const double x = xoffs + 0.5 * (2 * i - 1) * dx;
              const double w = ib_s_weight(
                  deltas,
                  pr,
                  (double [NDIMS]){px, py_, pz_},
                  (double [NDIMS]){x, y, z}
              );
              const double ux_p = pux + pvy * (z - pz_) - pvz * (y - py_);
              const double uy_p = puy + pvz * (x - px ) - pvx * (z - pz_);
              const double uz_p = puz + pvx * (y - py ) - pvy * (x - px );
              const double ux_f = + 0.5 * UX(i  , j  , k  ) + 0.5 * UX(i+1, j  , k  );
              const double uy_f = + 0.5 * UY(i  , j  , k  ) + 0.5 * UY(i  , j+1, k  );
              const double uz_f = + 0.5 * UZ(i  , j  , k  ) + 0.5 * UZ(i  , j  , k+1);
              const double fx = w * (ux_p - ux_f) / dt;
              const double fy = w * (uy_p - uy_f) / dt;
              const double fz = w * (uz_p - uz_f) / dt;
              IBFX(i, j, k) += fx * dt;
              IBFY(i, j, k) += fy * dt;
              IBFZ(i, j, k) += fz * dt;
              *fux -= + fx * cellsize / pm;
              *fuy -= + fy * cellsize / pm;
              *fuz -= + fz * cellsize / pm;
              *tvx -= - (z - pz_) * fy * cellsize / pmi;
              *tvx -= + (y - py ) * fz * cellsize / pmi;
              *tvy -= - (x - px ) * fz * cellsize / pmi;
              *tvy -= + (z - pz ) * fx * cellsize / pmi;
              *tvz -= - (y - py_) * fx * cellsize / pmi;
              *tvz -= + (x - px ) * fy * cellsize / pmi;
            }
          }
        }
      }
    }
#endif
    update_boundaries(domain, &ib->ibfx);
    update_boundaries(domain, &ib->ibfy);
#if NDIMS == 3
    update_boundaries(domain, &ib->ibfz);
#endif
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
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  double * ux = fluid->ux.data;
  double * uy = fluid->uy.data;
#if NDIMS == 3
  double * uz = fluid->uz.data;
#endif
  const double * ibfx = ib->ibfx.data;
  const double * ibfy = ib->ibfy.data;
#if NDIMS == 3
  const double * ibfz = ib->ibfz.data;
#endif
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      UX(i, j) +=
        + 0.5 * IBFX(i-1, j  )
        + 0.5 * IBFX(i  , j  );
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        UX(i, j, k) +=
          + 0.5 * IBFX(i-1, j  , k  )
          + 0.5 * IBFX(i  , j  , k  );
      }
    }
  }
#endif
  if(0 != fluid_update_boundaries_ux(domain, &fluid->ux)){
    return 1;
  }
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      UY(i, j) +=
        + 0.5 * IBFY(i  , j-1)
        + 0.5 * IBFY(i  , j  );
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        UY(i, j, k) +=
          + 0.5 * IBFY(i  , j-1, k  )
          + 0.5 * IBFY(i  , j  , k  );
      }
    }
  }
#endif
  if(0 != fluid_update_boundaries_uy(domain, &fluid->uy)){
    return 1;
  }
#if NDIMS == 3
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        UZ(i, j, k) +=
          + 0.5 * IBFZ(i  , j  , k-1)
          + 0.5 * IBFZ(i  , j  , k  );
      }
    }
  }
  if(0 != fluid_update_boundaries_uz(domain, &fluid->uz)){
    return 1;
  }
#endif
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
#if NDIMS == 2
  const size_t nvars = 3;
#else
  const size_t nvars = 6;
#endif
  double * buf = memory_calloc(nvars * nitems, sizeof(double));
  // pack
  for(size_t n = 0; n < nitems; n++){
    particle_t * p = particles + n;
#if NDIMS == 2
    buf[nvars * n + 0] = p->fux;
    buf[nvars * n + 1] = p->fuy;
    buf[nvars * n + 2] = p->tvz;
#else
    buf[nvars * n + 0] = p->fux;
    buf[nvars * n + 1] = p->fuy;
    buf[nvars * n + 2] = p->fuz;
    buf[nvars * n + 3] = p->tvx;
    buf[nvars * n + 4] = p->tvy;
    buf[nvars * n + 5] = p->tvz;
#endif
  }
  // sum up all
  MPI_Allreduce(MPI_IN_PLACE, buf, nvars * nitems, MPI_DOUBLE, MPI_SUM, comm_cart);
  // unpack
  for(size_t n = 0; n < nitems; n++){
    particle_t * p = particles + n;
#if NDIMS == 2
    p->fux = buf[nvars * n + 0];
    p->fuy = buf[nvars * n + 1];
    p->tvz = buf[nvars * n + 2];
#else
    p->fux = buf[nvars * n + 0];
    p->fuy = buf[nvars * n + 1];
    p->fuz = buf[nvars * n + 2];
    p->tvx = buf[nvars * n + 3];
    p->tvy = buf[nvars * n + 4];
    p->tvz = buf[nvars * n + 5];
#endif
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

