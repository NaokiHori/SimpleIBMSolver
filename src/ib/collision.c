#include <stdio.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include "sdecomp.h"
#include "memory.h"
#include "domain.h"
#include "ib.h"
#include "ib_solver.h"

static const double pi = 3.14159265358979324;

// Collision pair table
//
//    | 0  1    2     3 ...  N-2  N-1 <- n1
// ---+------------------------------
//  0 |    0    1     2 ...  N-3  N-2
//  1 |       N-1     N ... 2N-5 2N-4
//  2 |            2N-3 ... 3N-8 3N-7
// ...
//  ^ |
//  |
//  n0

static size_t get_l(
    const size_t nitems,
    const size_t n0
){
  // left -most index of the collision pair matrix for particle index "n0"
  return (2 * nitems - n0 - 1) * (n0    ) / 2    ;
}

static size_t get_r(
    const size_t nitems,
    const size_t n0
){
  // right-most index of the collision pair matrix for particle index "n0"
  return (2 * nitems - n0 - 2) * (n0 + 1) / 2 - 1;
}

static int get_particle_indices(
    const size_t nitems,
    const size_t n,
    size_t * n0,
    size_t * n1
){
  // convert index of collision pair "n" to the corresponding particle indices "n0" and "n1"
  for(*n0 = 0; ; (*n0)++){
    if(get_l(nitems, *n0) <= n && n <= get_r(nitems, *n0)){
      break;
    }
  }
  *n1 = n - get_l(nitems, *n0) + (*n0) + 1;
  return 0;
}

static int get_my_range(
    const domain_t * domain,
    const size_t nitems,
    size_t * n_min,
    size_t * n_max
){
  int nprocs_ = 0;
  int myrank_ = 0;
  sdecomp.get_comm_size(domain->info, &nprocs_);
  sdecomp.get_comm_rank(domain->info, &myrank_);
  size_t nprocs = nprocs_;
  size_t myrank = myrank_;
  *n_min = 0;
  for(size_t n = 0; n < myrank; n++){
    *n_min += (nitems + n) / nprocs;
  }
  *n_max = *n_min + (nitems + myrank) / nprocs;
  return 0;
}

static double harmonic_average(
    const double v0,
    const double v1
){
  const double retval = 1. / v0 + 1. / v1;
  return 2. / retval;
}

static int compute_collision_force_p_p(
    const domain_t * domain,
    const size_t index,
    particle_t * p0,
    particle_t * p1
){
  // correct periodicity
  double yoffset = 0.;
  {
    const double ly = domain->lengths[1];
    const double p0y = p0->y + p0->dy;
    const double p1y = p1->y + p1->dy;
    double minval = DBL_MAX;
    for(int yper = -1; yper <= 1; yper++){
      double val = fabs((p1y + ly * yper) - p0y);
      if(val < minval){
        minval = val;
        yoffset = ly * yper;
      }
    }
  }
#if NDIMS == 3
  double zoffset = 0.;
  {
    const double lz = domain->lengths[2];
    const double p0z = p0->z + p0->dz;
    const double p1z = p1->z + p1->dz;
    double minval = DBL_MAX;
    for(int zper = -1; zper <= 1; zper++){
      double val = fabs((p1z + lz * zper) - p0z);
      if(val < minval){
        minval = val;
        zoffset = lz * zper;
      }
    }
  }
#endif
  // check collision and compute collision force when needed
  const double r0 = p0->r;
  const double m0 = p0->m;
  const double x0 = p0->x + p0->dx;
  const double y0 = p0->y + p0->dy;
#if NDIMS == 3
  const double z0 = p0->z + p0->dz;
#endif
  const double r1 = p1->r;
  const double m1 = p1->m;
  const double x1 = p1->x + p1->dx;
  const double y1 = p1->y + p1->dy + yoffset;
#if NDIMS == 3
  const double z1 = p1->z + p1->dz + zoffset;
#endif
  // normal and overlap distance
  double nx = x1 - x0;
  double ny = y1 - y0;
#if NDIMS == 3
  double nz = z1 - z0;
#endif
  double norm = sqrt(
    + nx * nx
    + ny * ny
#if NDIMS == 3
    + nz * nz
#endif
  );
  // to avoid zero division just in case
  norm = fmax(norm, DBL_EPSILON);
  nx /= norm;
  ny /= norm;
#if NDIMS == 3
  nz /= norm;
#endif
  const double overlap_dist = r0 + r1 - norm;
  if(overlap_dist < 0.){
    // do not collide
    return 0;
  }
  // do collide
  // compute force
  // k: pre-factor (spring stiffness)
  const double m = harmonic_average(m0, m1);
  const double r = harmonic_average(r0, r1);
  const double k = m * (pi * pi) / (r * r);
  const double cfx = k * overlap_dist * nx;
  const double cfy = k * overlap_dist * ny;
#if NDIMS == 3
  const double cfz = k * overlap_dist * nz;
#endif
  p0->cfx[index] -= 1. / m0 * cfx;
  p1->cfx[index] += 1. / m1 * cfx;
  p0->cfy[index] -= 1. / m0 * cfy;
  p1->cfy[index] += 1. / m1 * cfy;
#if NDIMS == 3
  p0->cfz[index] -= 1. / m0 * cfz;
  p1->cfz[index] += 1. / m1 * cfz;
#endif
  return 0;
}

static int compute_collision_force_p_w(
    const double wallx,
    const size_t index,
    particle_t * p
){
  const double x = p->x + p->dx;
  const double r = p->r;
  const double m = p->m;
  double nx = wallx - x;
  const double norm = fabs(nx);
  nx /= norm;
  const double overlap_dist = r - norm;
  if(overlap_dist < 0.){
    return 0;
  }
  // compute force and torque
  // k: pre-factor (spring stiffness)
  const double k = m * (pi * pi) / (r * r);
  const double cfx = k * overlap_dist * nx;
  const double cfy = 0.;
#if NDIMS == 3
  const double cfz = 0.;
#endif
  p->cfx[index] -= 1. / m * cfx;
  p->cfy[index] -= 1. / m * cfy;
#if NDIMS == 3
  p->cfz[index] -= 1. / m * cfz;
#endif
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
#if NDIMS == 2
  const size_t nvars = 2;
#else
  const size_t nvars = 3;
#endif
  double * buf = memory_calloc(nvars * nitems, sizeof(double));
  // pack
  for(size_t n = 0; n < nitems; n++){
    particle_t * p = particles + n;
    buf[nvars * n + 0] = p->cfx[index];
    buf[nvars * n + 1] = p->cfy[index];
#if NDIMS == 3
    buf[nvars * n + 2] = p->cfz[index];
#endif
  }
  // sum up all
  MPI_Allreduce(MPI_IN_PLACE, buf, nvars * nitems, MPI_DOUBLE, MPI_SUM, comm_cart);
  // unpack
  for(size_t n = 0; n < nitems; n++){
    particle_t * p = particles + n;
    p->cfx[index] = buf[nvars * n + 0];
    p->cfy[index] = buf[nvars * n + 1];
#if NDIMS == 3
    p->cfz[index] = buf[nvars * n + 2];
#endif
  }
  memory_free(buf);
  return 0;
}

int ib_compute_collision_force(
    const domain_t * domain,
    const size_t index,
    ib_t * ib
){
  // NOTE: only the spring in the normal direction is considered for simplicity
  const double lx = domain->lengths[0];
  const size_t nps = ib->nitems;
  particle_t * particles = ib->particles;
  // reset collision force at this CN step
  for(size_t n = 0; n < nps; n++){
    particle_t * p = particles + n;
    p->cfx[index] = 0.;
    p->cfy[index] = 0.;
#if NDIMS == 3
    p->cfz[index] = 0.;
#endif
  }
  // particle-particle collisions
  {
    const size_t ncols = nps * (nps - 1) / 2;
    size_t n_min = 0;
    size_t n_max = 0;
    get_my_range(domain, ncols, &n_min, &n_max);
    for(size_t n = n_min; n < n_max; n++){
      size_t n0 = 0;
      size_t n1 = 0;
      get_particle_indices(nps, n, &n0, &n1);
      compute_collision_force_p_p(domain, index, particles + n0, particles + n1);
    }
  }
  // particle-wall collisions
  {
    size_t n_min = 0;
    size_t n_max = 0;
    get_my_range(domain, nps, &n_min, &n_max);
    for(size_t n = n_min; n < n_max; n++){
      compute_collision_force_p_w(0., index, particles + n);
      compute_collision_force_p_w(lx, index, particles + n);
    }
  }
  // synchronise computed forcings
  if(0 != synchronise_information(domain, index, ib)){
    return 1;
  }
  return 0;
}

