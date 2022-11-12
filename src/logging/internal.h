#if !defined(LOGGING_INTERNAL_H)
#define LOGGING_INTERNAL_H

// this header file should be used only inside src/logging/

#include "domain.h"
#include "fluid.h"
#include "temperature.h"

extern int check_divergence(const char fname[], const domain_t *domain, const double time, const fluid_t *fluid);
extern int check_momentum(const char fname[], const domain_t *domain, const double time, const fluid_t *fluid);
extern int check_energy(const char fname[], const domain_t *domain, const double time, const fluid_t *fluid, const temperature_t *temperature);
extern int check_nusselt(const char fname[], const domain_t *domain, const double time, const fluid_t *fluid, const temperature_t *temperature);

#endif // LOGGING_INTERNAL_H
