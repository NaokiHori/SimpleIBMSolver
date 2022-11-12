#if !defined(FLUID_CORRECT_VELOCITY_INTERNAL_H)
#define FLUID_CORRECT_VELOCITY_INTERNAL_H

#include "domain.h"
#include "fluid.h"

#if NDIMS == 2

extern int fluid_correct_velocity_ux(const domain_t * restrict domain, const double prefactor, fluid_t * restrict fluid);
extern int fluid_correct_velocity_uy(const domain_t * restrict domain, const double prefactor, fluid_t * restrict fluid);

#else // NDIMS == 3

extern int fluid_correct_velocity_ux(const domain_t * restrict domain, const double prefactor, fluid_t * restrict fluid);
extern int fluid_correct_velocity_uy(const domain_t * restrict domain, const double prefactor, fluid_t * restrict fluid);
extern int fluid_correct_velocity_uz(const domain_t * restrict domain, const double prefactor, fluid_t * restrict fluid);

#endif // NDIMS

#endif // FLUID_CORRECT_VELOCITY_INTERNAL_H
