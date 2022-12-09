#if !defined(FLUID_CORRECT_VELOCITY_INTERNAL_H)
#define FLUID_CORRECT_VELOCITY_INTERNAL_H

#include "domain.h"
#include "fluid.h"


extern int fluid_correct_velocity_ux(const domain_t * restrict domain, const double prefactor, fluid_t * restrict fluid);
extern int fluid_correct_velocity_uy(const domain_t * restrict domain, const double prefactor, fluid_t * restrict fluid);
extern int fluid_correct_velocity_uz(const domain_t * restrict domain, const double prefactor, fluid_t * restrict fluid);


#endif // FLUID_CORRECT_VELOCITY_INTERNAL_H
