#include "domain.h"
#include "fluid.h"
#include "arrays/fluid.h"


#if NDIMS == 2

/**
 * @brief update boundary values of x velocity
 * @param[in   ] domain : information about domain decomposition and size
 * @param[inout] ux     : x velocity
 * @return              : error code
 */
int fluid_update_boundaries_ux(const domain_t *domain, double *ux){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  /* ! update halo values ! 1 ! */
  domain_communicate_halo_ux_like(domain, ux);
  /* ! set boundary values ! 4 ! */
  for(int j = 1; j <= jsize; j++){
    UX(      1, j) = 0.; // impermeable
    UX(isize+1, j) = 0.; // impermeable
  }
  return 0;
}

#else // NDIMS == 3

/**
 * @brief update boundary values of x velocity
 * @param[in   ] domain : information about domain decomposition and size
 * @param[inout] ux     : x velocity
 * @return              : error code
 */
int fluid_update_boundaries_ux(const domain_t * restrict domain, double * restrict ux){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  /* ! update halo values ! 1 ! */
  domain_communicate_halo_ux_like(domain, ux);
  /* ! set boundary values ! 6 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      UX(      1, j, k) = 0.; // impermeable
      UX(isize+1, j, k) = 0.; // impermeable
    }
  }
  return 0;
}

#endif // NDIMS
