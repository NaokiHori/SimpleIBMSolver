#if !defined(INTEGRATE_H)
#define INTEGRATE_H

// main integrator
extern int integrate(
    const domain_t * domain,
    fluid_t * fluid,
    ib_t * ib,
    double * dt
);

#endif // INTEGRATE_H
