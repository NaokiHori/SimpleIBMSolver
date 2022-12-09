#include "domain.h"
#include "fluid.h"
#include "internal.h"



/**
 * @brief update velocity (can be non-solenoidal)
 * @param[in   ] domain : information related to MPI domain decomposition
 * @param[in   ] rkstep : Runge-Kutta step
 * @param[in   ] dt     : time step size
 * @param[inout] fluid  : velocity and pressure (in), velocity (out)
 * @return              : error code
 */
int fluid_update_velocity(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
  fluid_update_velocity_ux(domain, rkstep, dt, fluid);
  fluid_update_velocity_uy(domain, rkstep, dt, fluid);
  fluid_update_velocity_uz(domain, rkstep, dt, fluid);
  return 0;
}

/**
 * @brief destruct local buffers used in this source
 * @return : error code
 */
int fluid_update_velocity_finalise(void){
  fluid_update_velocity_finalise_ux();
  fluid_update_velocity_finalise_uy();
  fluid_update_velocity_finalise_uz();
  return 0;
}

