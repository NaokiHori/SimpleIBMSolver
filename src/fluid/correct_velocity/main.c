#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "arrays/fluid.h"
#include "internal.h"



/**
 * @brief correct non-solenoidal velocity using scalar potential \psi
 * @param[in   ] domain : information about domain decomposition and size
 * @param[in   ] rkstep : Runge-Kutta step
 * @param[in   ] dt     : time step size
 * @param[inout] fluid  : scalar potential \psi (in), velocity (out)
 * @return              : error code
 */
int fluid_correct_velocity(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
  const double gamma = RKCOEFS[rkstep].gamma;
  const double prefactor = gamma * dt;
  fluid_correct_velocity_ux(domain, prefactor, fluid);
  fluid_correct_velocity_uy(domain, prefactor, fluid);
  fluid_correct_velocity_uz(domain, prefactor, fluid);
  return 0;
}

