#if !defined(FLUID_H)
#define FLUID_H

#include "domain.h"
#include "structure.h"


/* ! definition of a structure fluid_t_ ! 15 !*/
/** @struct fluid_t_ (or fluid_t)
 *  @brief struct storing fluid-related variables
 *  @var ux, uy, uz             : velocity in each direction
 *  @var p, psi                 : pressure, scalar potential \psi
 *  @var srcuxa, srcuxb, srcuxg : Runge-Kutta source terms, ux
 *  @var srcuya, srcuyb, srcuyg : Runge-Kutta source terms, uy
 *  @var srcuza, srcuzb, srcuzg : Runge-Kutta source terms, uz
 */
struct fluid_t_ {
  double *ux, *uy, *uz;
  double *p, *psi;
  double *srcuxa, *srcuxb, *srcuxg;
  double *srcuya, *srcuyb, *srcuyg;
  double *srcuza, *srcuzb, *srcuzg;
};


// constructor and destructor
extern fluid_t *fluid_init(const domain_t * restrict domain);
extern int fluid_finalise(fluid_t * restrict fluid);

// compute right-hand-side of Runge-Kutta scheme
extern int fluid_compute_rhs(const domain_t * restrict domain, const int rkstep, fluid_t * restrict fluid, const temperature_t * restrict temperature);

// update velocity field and its clean-up function
extern int fluid_update_velocity(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid);
extern int fluid_update_velocity_finalise(void);

// compute scalar potential and its clean-up function
extern int fluid_compute_potential(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid);
extern int fluid_compute_potential_finalise(void);

// correct velocity field using scalar potential
extern int fluid_correct_velocity(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid);

// update pressure
extern int fluid_update_pressure(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid);

// exchange halos and impose boundary conditions
extern int fluid_update_boundaries_ux(const domain_t * restrict domain, double * restrict ux);
extern int fluid_update_boundaries_uy(const domain_t * restrict domain, double * restrict uy);
extern int fluid_update_boundaries_uz(const domain_t * restrict domain, double * restrict uz);
extern int fluid_update_boundaries_p (const domain_t * restrict domain, double * restrict  p);

#endif // FLUID_H
