#if !defined(FLUID_COMPUTE_RHS_INTERNAL)
#define FLUID_COMPUTE_RHS_INTERNAL

#include "domain.h"
#include "fluid.h"
#include "temperature.h"

extern int compute_rhs_ux(const domain_t * restrict domain, const int rkstep, fluid_t * restrict fluid, const temperature_t * restrict temperature);
extern int compute_rhs_uy(const domain_t * restrict domain, const int rkstep, fluid_t * restrict fluid);

#endif // FLUID_COMPUTE_RHS_INTERNAL
