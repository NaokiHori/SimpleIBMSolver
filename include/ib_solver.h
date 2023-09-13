#if !defined(IB_SOLVER_H)
#define IB_SOLVER_H

#include "domain.h"
#include "fluid.h"
#include "ib.h"

extern int ib_init(
    const char dirname_ic[],
    const domain_t * domain,
    ib_t * ib
);

extern int ib_reset_variables(
    ib_t * ib
);

extern int ib_exchange_momentum(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid,
    ib_t * ib
);

extern int ib_compute_inertia(
    const domain_t * domain,
    const size_t index,
    const fluid_t * fluid,
    ib_t * ib
);

extern int ib_compute_collision_force(
    const domain_t * domain,
    const size_t index,
    ib_t * ib
);

extern int ib_increment_particles(
    const size_t rkstep,
    double dt,
    ib_t * ib
);

extern int ib_update_particles(
    const domain_t * domain,
    ib_t * ib
);

extern int ib_save(
    const char dirname[],
    const domain_t * domain,
    const ib_t * ib
);

#endif // IB_SOLVER_H
