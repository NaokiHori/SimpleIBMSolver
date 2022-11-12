#if !defined(DOMAIN_INTERNAL_H)
#define DOMAIN_INTERNAL_H

#include "domain.h"

// coordinate.c
extern double *allocate_and_init_xf(const int isize, const double lx);
extern double *allocate_and_init_xc(const int isize, const double *xf);
extern double *allocate_and_init_dxf(const int isize, const double *xf);
extern double *allocate_and_init_dxc(const int isize, const double *xc);

// init.c
extern domain_t *domain_init(void);

// finalise.c
extern int domain_finalise(domain_t *domain);

// save.c
extern int domain_save(const char dirname[], const domain_t *domain);

#endif // DOMAIN_INTERNAL_H
