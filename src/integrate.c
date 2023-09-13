#include "runge_kutta.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "ib.h"
#include "ib_solver.h"
#include "integrate.h"
#include "decide_dt.h"

// integrate the equations for one time step
int integrate(
    const domain_t * domain,
    fluid_t * fluid,
    ib_t * ib,
    double * dt
){
  // decide time step size
  if(0 != decide_dt(domain, fluid, dt)){
    return 1;
  }
  // Runge-Kutta iterations
  // max iteration, should be three
  for(size_t rkstep = 0; rkstep < RKSTEPMAX; rkstep++){
    // reset ib-related right-hand-side variables
    if(0 != ib_reset_variables(ib)){
      return 1;
    }
    // predict flow field
    // compute right-hand-side terms of RK scheme
    if(0 != fluid_compute_rhs(domain, fluid)){
      return 1;
    }
    // compute advective terms inside particle at k step
    if(0 != ib_compute_inertia(domain, 0, fluid, ib)){
      return 1;
    }
    // compute collision force based on k-step position
    if(0 != ib_compute_collision_force(domain, 0, ib)){
      return 1;
    }
    // couple external factors: buoyancy force
    if(0 != fluid_couple_external_force(domain, fluid)){
      return 1;
    }
    // update flow field
    // NOTE: while the temperature is fully updated here,
    //   the velocity field is still non-solenoidal
    if(0 != fluid_predict_field(domain, rkstep, *dt, fluid)){
      return 1;
    }
    // compute hydrodynamic force between fluid and objects
    if(0 != ib_exchange_momentum(domain, rkstep, *dt, fluid, ib)){
      return 1;
    }
    // compute scalar potential
    // now the temperature field has been updated,
    //   while the velocity field is not divergence free
    //   and thus the following correction step is needed
    if(0 != fluid_compute_potential(domain, rkstep, *dt, fluid)){
      return 1;
    }
    // correct velocity field to satisfy mass conservation
    if(0 != fluid_correct_velocity(domain, rkstep, *dt, fluid)){
      return 1;
    }
    // update pressure
    if(0 != fluid_update_pressure(domain, rkstep, *dt, fluid)){
      return 1;
    }
    // update particles iteratively
    for(size_t substep = 0; substep < 4; substep++){
      // compute advective terms inside particle at k+1 step
      if(0 != ib_compute_inertia(domain, 1, fluid, ib)){
        return 1;
      }
      // compute collision force based on k-step position
      if(0 != ib_compute_collision_force(domain, 1, ib)){
        return 1;
      }
      // compute particle increments
      if(0 != ib_increment_particles(rkstep, *dt, ib)){
        return 1;
      }
    }
    // update particle velocity and position
    if(0 != ib_update_particles(domain, ib)){
      return 1;
    }
  }
  return 0;
}

