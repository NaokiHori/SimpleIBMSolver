#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
#include "config.h"
#include "param.h"
#include "array.h"
#include "sdecomp.h"
#include "domain.h"
#include "fluid.h"
#include "ib.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"

// overriden later using environment variables
static bool coefs_are_initialised = false;
static double coef_dt_adv = 0.;
static double coef_dt_dif = 0.;

/**
 * @brief decide time step size restricted by the advective terms
 * @param[in]  domain : information about domain decomposition and size
 * @param[in]  fluid  : velocity
 * @param[out] dt     : time step size
 * @return            : error code
 */
static int decide_dt_adv(
    const domain_t * domain,
    const fluid_t * fluid,
    double * restrict dt
){
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double dx = domain->dx;
  const double dy = domain->dy;
  const double * restrict ux = fluid->ux.data;
  const double * restrict uy = fluid->uy.data;
  // sufficiently small number to avoid zero division
  const double small = 1.e-8;
  // init with max possible dt
  *dt = 1.;
  // compute grid-size over velocity in x
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      double vel = fabs(UX(i, j)) + small;
      *dt = fmin(*dt, dx / vel);
    }
  }
  // compute grid-size over velocity in y
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      double vel = fabs(UY(i, j)) + small;
      *dt = fmin(*dt, dy / vel);
    }
  }
  // compute grid-size over velocity in z
  // unify result, multiply safety factor
  MPI_Allreduce(MPI_IN_PLACE, dt, 1, MPI_DOUBLE, MPI_MIN, comm_cart);
  *dt *= coef_dt_adv;
  return 0;
}

/**
 * @brief decide time step size restricted by the diffusive terms
 * @param[in]  domain      : grid size
 * @param[in]  diffusivity : fluid / temperature diffusivity
 * @param[out] dt          : time step size
 * @return                 : error code
 */
static int decide_dt_dif(
    const domain_t * domain,
    const double diffusivity,
    double * dt
){
  const double dx = domain->dx;
  const double dy = domain->dy;
  double grid_sizes[NDIMS] = {0.};
  grid_sizes[0] = dx;
  grid_sizes[1] = dy;
  // compute diffusive constraints
  for(int dim = 0; dim < NDIMS; dim++){
    dt[dim] = coef_dt_dif / diffusivity * 0.5 / NDIMS * pow(grid_sizes[dim], 2.);
  }
  return 0;
}

static int decide_dt_ib(
    const domain_t * domain,
    const ib_t * ib,
    double * dt
){
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  const double dx = domain->dx;
  const double dy = domain->dy;
  const size_t nitems = ib->nitems;
  const particle_t * particles = ib->particles;
  // sufficiently small number to avoid zero division
  const double small = 1.e-8;
  // init with max possible dt
  *dt = 1.;
  // compute grid-size over velocity in x
  for(size_t n = 0; n < nitems; n++){
    const particle_t * p = particles + n;
    double ux = fabs(p->ux) + p->r * fabs(p->vz) + small;
    double uy = fabs(p->uy) + p->r * fabs(p->vz) + small;
    *dt = fmin(*dt, dx / ux);
    *dt = fmin(*dt, dy / uy);
  }
  return 0;
}

/**
 * @brief decide time step size which can integrate the equations stably
 * @param[in]  domain : information about domain decomposition and size
 * @param[in]  fluid  : velocity and diffusivities
 * @param[in]  ib     : particle velocities
 * @param[out]        : time step size
 * @return            : (success) 0
 *                    : (failure) non-zero value
 */
int decide_dt(
    const domain_t * domain,
    const fluid_t * fluid,
    const ib_t * ib,
    double * dt
){
  if(!coefs_are_initialised){
    if(0 != config.get_double("coef_dt_adv", &coef_dt_adv)) return 1;
    if(0 != config.get_double("coef_dt_dif", &coef_dt_dif)) return 1;
    coefs_are_initialised = true;
    const int root = 0;
    int myrank = root;
    sdecomp.get_comm_rank(domain->info, &myrank);
    if(root == myrank){
      printf("coefs: (adv) % .3e, (dif) % .3e\n", coef_dt_adv, coef_dt_dif);
    }
  }
  // compute advective and diffusive constraints
  double dt_adv[1] = {0.};
  double dt_dif_m[NDIMS] = {0.};
  double dt_dif_t[NDIMS] = {0.};
  double dt_ib[1] = {0.};
  decide_dt_adv(domain, fluid,        dt_adv  );
  decide_dt_dif(domain, fluid->m_dif, dt_dif_m);
  decide_dt_dif(domain, fluid->t_dif, dt_dif_t);
  decide_dt_ib (domain, ib, dt_ib);
  // choose smallest value as dt
  // advection
  *dt = dt_adv[0];
  // diffusion, momentum
  if(!param_m_implicit_x){
    *dt = fmin(*dt, dt_dif_m[0]);
  }
  if(!param_m_implicit_y){
    *dt = fmin(*dt, dt_dif_m[1]);
  }
  // diffusion, temperature
  if(!param_t_implicit_x){
    *dt = fmin(*dt, dt_dif_t[0]);
  }
  if(!param_t_implicit_y){
    *dt = fmin(*dt, dt_dif_t[1]);
  }
  // ib-related
  *dt = fmin(*dt, dt_ib[0]);
  return 0;
}

