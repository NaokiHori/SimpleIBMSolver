#include <stdio.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include "sdecomp.h"
#include "common.h"
#include "domain.h"
#include "ib.h"


/*
 * Collision pair matrix
 *
 *    | 0  1    2     3 ...  N-2  N-1  <- n1
 * ---+------------------------------
 *  0 |    0    1     2 ...  N-3  N-2
 *  1 |       N-1     N ... 2N-5 2N-4
 *  2 |            2N-3 ... 3N-8 3N-7
 * ...
 *  ^ |
 *  |
 *  n0
 */

static int get_l(const int n_particles, const int n0){
  // left-most index of the collision pair matrix for particle index "n0"
  return (2*n_particles-n0-1)*(n0  )/2  ;
}

static int get_r(const int n_particles, const int n0){
  // right-most index of the collision pair matrix for particle index "n0"
  return (2*n_particles-n0-2)*(n0+1)/2-1;
}

static int get_particle_indices(const int n_particles, const int n, int *n0, int *n1){
  // convert index of collision pair "n" to the corresponding particle indices "n0" and "n1"
  for(*n0 = 0; ; (*n0)++){
    if(get_l(n_particles, *n0) <= n && n <= get_r(n_particles, *n0)){
      break;
    }
  }
  *n1 = n-get_l(n_particles, *n0)+(*n0)+1;
  return 0;
}

static int get_my_range(const domain_t *domain, const int n_total, int *n_min, int *n_max){
  const MPI_Comm comm = domain->sdecomp->comm_cart;
  int nprocs, myrank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myrank);
  *n_min = sdecomp_kernel_get_offset(n_total, nprocs, myrank);
  *n_max = sdecomp_kernel_get_mysize(n_total, nprocs, myrank) + *n_min;
  return 0;
}

static double harmonic_average(const double v0, const double v1){
  double retval = 0.5*(1./v0+1./v1);
  return 1./retval;
}


typedef struct circle_t_ {
  // center
  double x, y;
  // radius
  double r;
} circle_t;

static int compute_collision_force_p_p(const domain_t *domain, const int cnstep, particle_t *p0, particle_t *p1){
  /* correct periodicity */
  double yoffset = 0.;
  {
    const double ly = domain->lengths[1];
    const double p0y = p0->y + p0->dy;
    const double p1y = p1->y + p1->dy;
    double minval = DBL_MAX;
    for(int periodic = -1; periodic <= 1; periodic++){
      double val = fabs((p1y + ly * periodic) - p0y);
      if(val < minval){
        minval = val;
        yoffset = ly * periodic;
      }
    }
  }
  /* check collision and compute collision force when needed */
  const double x0 = p0->x + p0->dx;
  const double y0 = p0->y + p0->dy;
  const double r0 = p0->r;
  const double x1 = p1->x + p1->dx;
  const double y1 = p1->y + p1->dy + yoffset;
  const double r1 = p1->r;
  const double p0mass = ib_compute_mass(1., p0->r);
  const double p1mass = ib_compute_mass(1., p1->r);
  // normal and overlap distance
  double nx = x1 - x0;
  double ny = y1 - y0;
  double norm = sqrt(
    + pow(nx, 2.)
    + pow(ny, 2.)
  );
  // to avoid zero division just in case
  norm = fmax(norm, DBL_EPSILON);
  nx /= norm;
  ny /= norm;
  double overlap_dist = r0 + r1 - norm;
  if(overlap_dist < 0.){
    // do not collide
    return 0;
  }
  // do collide
  // compute force
  // k: pre-factor (spring stiffness)
  double k;
  {
    const double m = harmonic_average(p0mass, p1mass);
    const double r = harmonic_average(    r0,     r1);
    k = m * pow(M_PI, 2.)/pow(r, 2.);
  }
  double cfx = k * overlap_dist * nx;
  double cfy = k * overlap_dist * ny;
  // force in x
  p0->cfx[cnstep] += 1./p0mass*(-cfx);
  p1->cfx[cnstep] += 1./p1mass*(+cfx);
  // force in y
  p0->cfy[cnstep] += 1./p0mass*(-cfy);
  p1->cfy[cnstep] += 1./p1mass*(+cfy);
  return 0;
}

static int compute_collision_force_p_w(const double wallx, const int cnstep, particle_t *p){
  /* check collision and compute collision force when needed */
  const double x = p->x + p->dx;
  const double r = p->r;
  double nx = wallx-x;
  double norm = fabs(nx);
  nx /= norm;
  double overlap_dist = r - norm;
  if(overlap_dist < 0.){
    return 0;
  }
  // compute force and torque
  const double m = ib_compute_mass(1., p->r);
  // k: pre-factor (spring stiffness)
  double k = m * pow(M_PI, 2.)/pow(r, 2.);
  double cfx = k * overlap_dist * nx;
  double cfy = 0.;
  // force in x
  p->cfx[cnstep] += 1./m*(-cfx);
  // force in y
  p->cfy[cnstep] += 1./m*(-cfy);
  return 0;
}

int ib_compute_collision_force(const domain_t *domain, const int cnstep, ib_t *ib){
  /*
   * NOTE: only the spring in the normal direction is considered for simplicity
   * Although this is sufficient to avoid over-penetrations between particles,
   *   obviously not collect from a physical perspective
   */
  const double lx = domain->lengths[0];
  const int n_particles = ib->n_particles;
  particle_t **particles = ib->particles;
  // reset collision force at this CN step
  for(int n = 0; n < n_particles; n++){
    particle_t *p = particles[n];
    p->cfx[cnstep] = 0.;
    p->cfy[cnstep] = 0.;
  }
  // particle-particle collisions
  {
    const int n_total = n_particles*(n_particles-1)/2;
    int n_min, n_max;
    get_my_range(domain, n_total, &n_min, &n_max);
    for(int n = n_min; n < n_max; n++){
      int n0, n1;
      get_particle_indices(n_particles, n, &n0, &n1);
      compute_collision_force_p_p(domain, cnstep, particles[n0], particles[n1]);
    }
  }
  // particle-wall collisions
  {
    double wall_locations[2] = {0., lx};
    for(int wall_index = 0; wall_index < 2; wall_index++){
      double wall_location = wall_locations[wall_index];
      int n_min, n_max;
      get_my_range(domain, n_particles, &n_min, &n_max);
      for(int n = n_min; n < n_max; n++){
        compute_collision_force_p_w(wall_location, cnstep, particles[n]);
      }
    }
  }
  // synchronise computed forcings
  {
    // prepare message buffer
    double *buf = common_calloc(2*n_particles, sizeof(double));
    // pack
    for(int n = 0; n < n_particles; n++){
      particle_t *p = particles[n];
      buf[2*n+0] = p->cfx[cnstep];
      buf[2*n+1] = p->cfy[cnstep];
    }
    // sum up all
    MPI_Allreduce(MPI_IN_PLACE, buf, 2*n_particles, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // unpack
    for(int n = 0; n < n_particles; n++){
      particle_t *p = particles[n];
      p->cfx[cnstep] = buf[2*n+0];
      p->cfy[cnstep] = buf[2*n+1];
    }
    common_free(buf);
  }
  return 0;
}

