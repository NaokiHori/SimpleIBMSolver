#if !defined(SAVE_H)
#define SAVE_H

#include "domain.h"
#include "fluid.h"
#include "temperature.h"
#include "ib.h"

// next time to trigger output
extern double save_next;

extern int save(const domain_t *domain, const int step, const double time, const fluid_t *fluid, const temperature_t *temperature, const ib_t *ib);

#endif // SAVE_H
