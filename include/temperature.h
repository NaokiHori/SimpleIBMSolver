#if !defined(TEMPERATURE_H)
#define TEMPERATURE_H

#include "domain.h"
#include "structure.h"

/* ! definition of a structure temperature_t_ ! 11 ! */
/** @struct temperature_t_ (or temperature_t)
 *  @brief struct storing temperature-related variables
 *  @var temp                         : temperature
 *  @var tempforcex                   : buoyancy force
 *  @var srctempa, srctempb, srctempg : Runge-Kutta source terms, temp
 */
struct temperature_t_ {
  double *temp;
  double *tempforcex;
  double *srctempa, *srctempb, *srctempg;
};

// constructor and destructor
extern temperature_t *temperature_init(const domain_t * restrict domain);
extern int temperature_finalise(temperature_t * restrict temperature);

// compute buoyancy force
extern int temperature_compute_force(const domain_t * restrict domain, temperature_t * restrict temperature);

// compute right-hand-side of Runge-Kutta scheme
extern int temperature_compute_rhs(const domain_t * restrict domain, const int rkstep, const fluid_t * restrict fluid, temperature_t * restrict temperature);

// update temperature field and its clean-up function
extern int temperature_update_temp(const domain_t * restrict domain, const int rkstep, const double dt, temperature_t * restrict temperature);
extern int temperature_update_temp_finalise(void);

// exchange halos and impose boundary conditions
extern int temperature_update_boundaries_temp(const domain_t * restrict domain, double * restrict temp);

#endif // TEMPERATURE_H
